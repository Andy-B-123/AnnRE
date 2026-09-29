"""Portable prepare -> extract (one library/task) -> report workflow."""
import argparse
import csv
import gzip
import hashlib
import html
import json
import os
import platform
from collections import Counter, defaultdict
from pathlib import Path
import pysam
from . import __version__
from .core import (load_gtf, ident, union, geometry, overlap, observation, rejection,
                   motif, CANONICAL, classify_reference, preserved)


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, rows, fields):
    with open(path, 'w') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def dump(path, obj):
    temp = Path(str(path)+'.tmp')
    temp.write_text(json.dumps(obj, indent=2)+'\n')
    temp.replace(path)


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1024*1024), b''):
            h.update(b)
    return h.hexdigest()


def code_hashes():
    return {p.name: digest(p) for p in Path(__file__).parent.glob('*.py')}


def prepare(args):
    cfg = json.loads(Path(args.config).read_text())
    base = Path(args.config).resolve().parent
    for key in ('genome', 'reference', 'samples'):
        cfg[key] = str((base / cfg[key]).resolve())
    defaults = dict(min_mapq=20, max_clip_fraction=.1, max_edit_fraction=.05,
                    min_reads=3, min_samples=2, min_flank_reads=10, anchor=20,
                    flank_bases=50, observation_mode='pacbio-zmw', strand_mode='alignment',
                    context_bp=5000)
    for k, v in defaults.items():
        cfg.setdefault(k, v)
    if cfg['strand_mode'] not in ('alignment', 'ts') or cfg['observation_mode'] not in ('pacbio-zmw', 'read-name'):
        raise ValueError('Invalid strand/observation mode')
    if min(cfg[k] for k in ('anchor', 'min_reads', 'min_samples', 'min_flank_reads', 'flank_bases')) < 1:
        raise ValueError('Thresholds must be positive')
    if not 0 <= cfg['min_mapq'] < 255 or not all(0 <= cfg[k] <= 1 for k in ('max_clip_fraction', 'max_edit_fraction')):
        raise ValueError('Invalid quality thresholds')
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if (out/'prepared.json').exists():
        raise ValueError('Output already prepared; choose a new run directory')
    all_models = load_gtf(cfg['reference'])
    requested = set(cfg.get('genes', []))
    available = {m['gene_id'] for m in all_models.values()}
    if requested-available:
        raise ValueError(f'Unknown requested genes: {sorted(requested-available)}')
    selected = [m for m in all_models.values() if not requested or m['gene_id'] in requested]
    windows = defaultdict(list)
    with pysam.FastaFile(cfg['genome']) as fa:
        lengths = dict(zip(fa.references, fa.lengths))
        for m in all_models.values():
            if m['chrom'] not in lengths or m['exons'][-1][1] > lengths[m['chrom']]:
                raise ValueError('Reference annotation is incompatible with FASTA dictionary')
        for m in selected:
            windows[m['chrom']].append((max(0, m['exons'][0][0]-cfg['context_bp']),
                                       min(lengths[m['chrom']], m['exons'][-1][1]+cfg['context_bp'])))
    windows = {k: union(v) for k, v in windows.items()}
    # Whole neighboring models are retained as context, but reads are fetched only in selected windows.
    models = [m for m in all_models.values() if any(m['exons'][0][0] < e and m['exons'][-1][1] > s for s, e in windows.get(m['chrom'], []))]
    refs = {}
    for m in selected:
        for left, right in zip(m['exons'], m['exons'][1:]):
            key = (m['gene_id'], m['chrom'], left[1], right[0], m['strand'])
            eid = ident(key)
            r = refs.setdefault(eid, dict(event_id=eid, kind='reference_junction', genes=m['gene_id'],
                 chrom=m['chrom'], start=left[1], end=right[0], strand=m['strand'], transcripts=[], flank_pairs=[]))
            r['transcripts'].append(m['transcript_id'])
            r['flank_pairs'].append((left, right))
    samples = read_tsv(cfg['samples'])
    if not samples or any(not all(r.get(k) for k in ('library_id', 'sample_id', 'bam')) for r in samples):
        raise ValueError('Sample TSV requires nonempty library_id, sample_id, bam')
    if len({r['library_id'] for r in samples}) != len(samples):
        raise ValueError('Duplicate library IDs')
    by_sample = {}
    for row in samples:
        row['bam'] = str((Path(cfg['samples']).parent/row['bam']).resolve())
        info = tuple(row.get(k, '') for k in ('tissue', 'stage', 'instar'))
        if row['sample_id'] in by_sample and by_sample[row['sample_id']] != info:
            raise ValueError('Conflicting metadata within sample ID')
        by_sample[row['sample_id']] = info
    alternates = {}
    alternate_input_qc = {}
    for label, path in cfg.get('alternates', {}).items():
        path = str((base/path).resolve())
        cfg['alternates'][label] = path
        alternate_input_qc[label] = {}
        alternates[label] = [m for m in load_gtf(path, allow_unstranded=True, stats=alternate_input_qc[label]).values() if any(m['exons'][0][0] < e and m['exons'][-1][1] > s for s, e in windows.get(m['chrom'], []))]
    sources = {k: {'path': cfg[k], 'sha256': digest(cfg[k])} for k in ('reference', 'samples')}
    sources['genome'] = {'path': cfg['genome'], 'size': os.path.getsize(cfg['genome']),
                         'fai_sha256': digest(cfg['genome']+'.fai'),
                         'identity_limit': 'Dictionary/FAI validation; not a genome sequence checksum'}
    for label, path in cfg.get('alternates', {}).items():
        sources['alternate:'+label] = {'path': path, 'sha256': digest(path)}
    data = dict(version=__version__, code_sha256=code_hashes(), config=cfg, windows=windows, models=models,
                reference_events=list(refs.values()), alternates=alternates, alternate_input_qc=alternate_input_qc, samples=samples,
                contigs=lengths, provenance=sources, python=platform.python_version(), pysam=pysam.__version__)
    dump(out/'prepared.json', data)
    print(f'Prepared {len(refs)} reference junctions, {len(selected)} selected models, {len(samples)} libraries')


def extract(args):
    out = Path(args.run)
    prep_path = out/'prepared.json'
    p = json.loads(prep_path.read_text())
    cfg = p['config']
    if p['code_sha256'] != code_hashes():
        raise ValueError('Code changed since prepare; use a new run directory')
    row = p['samples'][args.library]
    dest = out/'evidence'/f'{args.library:04d}'
    dest.mkdir(parents=True, exist_ok=True)
    if (dest/'done.json').exists():
        raise ValueError('Library already complete; remove its evidence directory to rerun')
    genes = defaultdict(list)
    by_chr = defaultdict(list)
    for m in p['models']:
        genes[(m['chrom'], m['strand'], m['gene_id'])].extend(m['exons'])
        by_chr[(m['chrom'], m['strand'])].append(m)
    genes = {k: union(v) for k, v in genes.items()}
    ref_by_chr = defaultdict(list)
    for r in p['reference_events']:
        gx = genes[(r['chrom'], r['strand'], r['genes'])]
        r['left_context'] = [(a, min(b, r['start'])) for a, b in gx if a < r['start']]
        r['right_context'] = [(max(a, r['end']), b) for a, b in gx if b > r['end']]
        ref_by_chr[(r['chrom'], r['strand'])].append(r)
    events = {r['event_id']: r for r in p['reference_events']}
    qc = Counter()
    chains = {}
    records = dest/'observations.jsonl.gz'
    tmp = Path(str(records)+'.tmp')
    anchor = cfg['anchor']
    with pysam.AlignmentFile(row['bam'], 'rb') as bam, gzip.open(tmp, 'wt') as handle:
        if not bam.has_index():
            raise ValueError('BAM must be indexed')
        bam_lengths = dict(zip(bam.references, bam.lengths))
        if any(bam_lengths.get(c) != n for c, n in p['contigs'].items()):
            raise ValueError('BAM/FASTA contig dictionary mismatch')
        def emit(eid, label, obs, read, detail=''):
            cid = ident((chrom, strand, gaps))
            chains.setdefault(cid, dict(chrom=chrom, strand=strand, gaps=gaps, exons=exons, example_read=read.query_name))
            handle.write(json.dumps([eid, label, obs, read.query_name, detail, cid], separators=(',', ':'))+'\n')
        for chrom, windows in p['windows'].items():
            previous_end = -1
            for start, end in windows:
                for read in bam.fetch(chrom, start, end):
                    if read.reference_start < previous_end:
                        continue
                    qc['records_examined'] += 1
                    bad = rejection(read, cfg)
                    if bad:
                        qc.update('rejected:'+b for b in bad)
                        continue
                    try:
                        obs = observation(read.query_name, cfg['observation_mode'])
                    except ValueError:
                        qc['rejected:unresolved_observation'] += 1
                        continue
                    strand = '-' if read.is_reverse else '+'
                    if cfg['strand_mode'] == 'ts':
                        if not read.has_tag('ts') or read.get_tag('ts') not in ('+', '-'):
                            qc['rejected:missing_ts'] += 1
                            continue
                        # minimap2 ts is transcript strand relative to aligned strand.
                        if read.get_tag('ts') == '-':
                            strand = '+' if strand == '-' else '-'
                    qc['credible_records'] += 1
                    blocks, gaps, exons = geometry(read.reference_start, read.cigartuples)
                    strict = {tuple(g) for g in gaps if overlap(blocks, (g[0]-anchor, g[0])) == anchor and overlap(blocks, (g[1], g[1]+anchor)) == anchor}
                    for r in ref_by_chr[(chrom, strand)]:
                        if read.reference_end <= min(x[0] for x in r['left_context']) or read.reference_start >= max(x[1] for x in r['right_context']):
                            continue
                        flags = set()
                        for l, rr in r['flank_pairs']:
                            lok = overlap(blocks, l) >= min(cfg['flank_bases'], l[1]-l[0])
                            rok = overlap(blocks, rr) >= min(cfg['flank_bases'], rr[1]-rr[0])
                            if lok: flags.add('left')
                            if rok: flags.add('right')
                            if lok and rok: flags.add('bridge')
                        if (r['start'], r['end']) in strict:
                            flags.add('exact')
                        context_bridge = all(sum(overlap(blocks, x) for x in r[k]) >= cfg['flank_bases'] for k in ('left_context', 'right_context'))
                        if context_bridge:
                            flags.add('gene_context_bridge')
                        for flag in flags:
                            emit(r['event_id'], flag, obs, read)
                        for s, e in sorted(strict):
                            if (s, e) != (r['start'], r['end']) and (s == r['start'] or e == r['end']):
                                emit(r['event_id'], 'shared_boundary_alternative', obs, read, f'{s}:{e}')
                        if 'bridge' in flags or context_bridge:
                            for s, e in sorted(strict):
                                if (s, e) != (r['start'], r['end']) and s < r['end'] and e > r['start']:
                                    emit(r['event_id'], 'alternative', obs, read, f'{s}:{e}')
                    active = {gid: ex for (c, st, gid), ex in genes.items() if c == chrom and st == strand and any(overlap(blocks, x) for x in ex)}
                    def owners(segment):
                        return {g for g, xx in active.items() if sum(overlap([segment], x) for x in xx) >= anchor}
                    model_hits = defaultdict(list)
                    for m in by_chr[(chrom, strand)]:
                        if preserved(m, gaps, blocks, anchor):
                            model_hits[m['gene_id']].append(m['transcript_id'])
                    # Novel exon identities include both exact flanking introns.
                    for j in range(1, len(exons)-1):
                        exon = exons[j]
                        a, b = tuple(gaps[j-1]), tuple(gaps[j])
                        if a not in strict or b not in strict:
                            continue
                        common = owners(exons[j-1]) & owners(exons[j+1])
                        for gene in sorted(common):
                            if any(tuple(exon) == tuple(x) for m in by_chr[(chrom, strand)] if m['gene_id'] == gene for x in m['exons']):
                                continue
                            key = ('novel_exon_context', gene, chrom, strand, a, exon, b)
                            eid = ident(key)
                            events[eid] = dict(event_id=eid, kind=key[0], genes=gene, chrom=chrom, start=exon[0], end=exon[1], strand=strand, context=[a, b])
                            emit(eid, 'support', obs, read)
                    for j, (s, e) in enumerate(gaps):
                        if (s, e) not in strict:
                            continue
                        for left in sorted(owners(exons[j])):
                            for right in sorted(owners(exons[j+1])):
                                if left == right:
                                    continue
                                if any(sum(overlap(blocks, x) for x in active[g]) < cfg['flank_bases'] for g in (left, right)):
                                    continue
                                eid = ident(('cross_gene_junction', left, right, chrom, s, e, strand))
                                events[eid] = dict(event_id=eid, kind='cross_gene_junction', genes=left+';'+right, chrom=chrom, start=s, end=e, strand=strand)
                                emit(eid, 'support', obs, read)
                                if model_hits[left] and model_hits[right]:
                                    emit(eid, 'preserve_both', obs, read, ';'.join(sorted(model_hits[left]+model_hits[right])))
                previous_end = end
    tmp.replace(records)
    dump(dest/'events.json', events)
    dump(dest/'chains.json', chains)
    dump(dest/'done.json', dict(library_index=args.library, library_id=row['library_id'], sample_id=row['sample_id'],
                               prepared_sha256=digest(prep_path), bam=row['bam'], bam_size=os.path.getsize(row['bam']),
                               bam_mtime_ns=os.stat(row['bam']).st_mtime_ns, qc=dict(qc),
                               observations_sha256=digest(records), events_sha256=digest(dest/'events.json'), chains_sha256=digest(dest/'chains.json')))
    print(f'{row["library_id"]}: {dict(qc)}')


def model_matches(event, model):
    if (event['chrom'], event['strand']) != (model['chrom'], model['strand']):
        return False
    gaps = [tuple(x) for x in model['chain']]
    if event['kind'] != 'novel_exon_context':
        return (event['start'], event['end']) in gaps
    a, b = [tuple(x) for x in event['context']]
    return any(gaps[i:i+2] == [a, b] for i in range(len(gaps)-1))


def plot_event(event, models, alternates, chains, counts):
    """Native SVG, genomic coordinates on one linear scale; evidence, not an IGV image."""
    tracks = []
    selected_genes = event['genes'].split(';')
    for m in models:
        if m['gene_id'] in selected_genes:
            tracks.append(('Ref '+m['transcript_id'], m['exons'], '#506d86'))
    if tracks:
        lo = min(x[0] for _, ex, _ in tracks for x in ex)
        hi = max(x[1] for _, ex, _ in tracks for x in ex)
    else:
        lo, hi = event['start']-1000, event['end']+1000
    alt_omitted = 0
    for source, mm in alternates.items():
        matching = [m for m in mm if m['chrom'] == event['chrom'] and m['strand'] == event['strand'] and m['exons'][0][0] < hi and m['exons'][-1][1] > lo]
        alt_omitted += max(0, len(matching)-6)
        for m in matching[:6]:
            tracks.append((source+' '+m['transcript_id'], m['exons'], '#8f6daf'))
    for cid, n in counts.most_common(3):
        tracks.append((f'Observed chain: {n} supporting reads', chains[cid]['exons'], '#159079'))
    # Extent is enlarged to avoid clipping displayed models.
    lo = min([event['start']]+[a for _, ex, _ in tracks for a, b in ex])
    hi = max([event['end']]+[b for _, ex, _ in tracks for a, b in ex])
    x = lambda p: 310+860*(p-lo)/max(1, hi-lo)
    h = 110+30*len(tracks)
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1200 {h}" role="img" aria-label="Annotation and observed splice chains">',
           '<rect width="100%" height="100%" fill="white"/>',
           f'<text x="12" y="22" font-family="sans-serif" font-size="14">{html.escape(event["genes"])} · {html.escape(event["strand"])} strand</text>',
           f'<text x="310" y="43" font-size="12">{html.escape(event["chrom"])}:{lo+1:,}–{hi:,} (1-based display)</text>',
           f'<rect x="{x(event["start"]):.2f}" y="50" width="{max(2,x(event["end"])-x(event["start"])):.2f}" height="{h-90}" fill="#fde4bf"/>']
    for i, (label, exons, color) in enumerate(tracks):
        y = 68+30*i
        svg.append(f'<text x="12" y="{y+4}" font-size="11" font-family="sans-serif">{html.escape(label)}</text>')
        svg.append(f'<line x1="{x(exons[0][0]):.2f}" x2="{x(exons[-1][1]):.2f}" y1="{y}" y2="{y}" stroke="{color}"/>')
        for a, b in exons:
            svg.append(f'<rect x="{x(a):.2f}" y="{y-5}" width="{max(1,x(b)-x(a)):.2f}" height="10" fill="{color}"/>')
    svg.append(f'<text x="12" y="{h-12}" font-size="11">Orange: candidate interval. Observed chains use representative ends; abundance is event-associated. Alternate models omitted: {alt_omitted}.</text></svg>')
    return ''.join(svg)


def report(args):
    out = Path(args.run)
    p = json.loads((out/'prepared.json').read_text())
    cfg = p['config']
    if p['code_sha256'] != code_hashes():
        raise ValueError('Code changed since prepare; use a new run directory')
    phash = digest(out/'prepared.json')
    all_events = {r['event_id']: r for r in p['reference_events']}
    support = defaultdict(set)
    raw = Counter()
    chain_support = defaultdict(set)
    chains = {}
    manifests = []
    # A report is never silently produced from an incomplete array.
    for i, row in enumerate(p['samples']):
        d = out/'evidence'/f'{i:04d}'
        if not (d/'done.json').exists():
            raise ValueError(f'Incomplete run: missing library {i} ({row["library_id"]})')
        done = json.loads((d/'done.json').read_text())
        if done['prepared_sha256'] != phash or done['library_id'] != row['library_id']:
            raise ValueError('Evidence provenance mismatch')
        for name, key in [('observations.jsonl.gz', 'observations_sha256'), ('events.json', 'events_sha256'), ('chains.json', 'chains_sha256')]:
            if digest(d/name) != done[key]:
                raise ValueError(f'Evidence checksum mismatch: {d/name}')
        manifests.append(done)
        all_events.update(json.loads((d/'events.json').read_text()))
        chains.update(json.loads((d/'chains.json').read_text()))
        with gzip.open(d/'observations.jsonl.gz', 'rt') as f:
            for line in f:
                eid, label, obs, read_id, detail, cid = json.loads(line)
                key = (eid, label, detail if label in ('alternative', 'shared_boundary_alternative') else '')
                support[key].add((row['sample_id'], obs))
                raw[key] += 1
                if label in ('support', 'bridge', 'gene_context_bridge', 'left', 'right', 'exact'):
                    chain_support[(eid, cid)].add((row['sample_id'], obs))
    count = lambda s: (len(s), len({sample for sample, obs in s}))
    rows, sample_rows, alt_rows = [], [], []
    alternative_sets = defaultdict(dict)
    shared_boundary_sets = defaultdict(dict)
    for (eid, label, detail), ss in support.items():
        if label == 'alternative':
            alternative_sets[eid][detail] = ss
        elif label == 'shared_boundary_alternative':
            shared_boundary_sets[eid][detail] = ss
    meta = {r['sample_id']: r for r in p['samples']}
    with pysam.FastaFile(cfg['genome']) as fa:
        for eid, event in sorted(all_events.items()):
            r = {k: event[k] for k in ('event_id', 'kind', 'genes', 'chrom', 'start', 'end', 'strand')}
            sets = {k: support[(eid, k, '')] for k in ('exact', 'bridge', 'gene_context_bridge', 'left', 'right', 'support', 'preserve_both')}
            alt = alternative_sets[eid]
            shared = shared_boundary_sets[eid]
            if event['kind'] == 'reference_junction':
                coexp = {s for s, o in sets['left']} & {s for s, o in sets['right']}
                category = classify_reference(count(sets['exact']), [count(s) for s in alt.values()], count(sets['bridge'] | sets['gene_context_bridge']), len(sets['left']), len(sets['right']), len(coexp), cfg)
                if category not in ('EXACT_JUNCTION_SUPPORTED', 'ALTERNATIVE_CONNECTION_SUPPORTED') and any(count(ss)[0] >= cfg['min_reads'] and count(ss)[1] >= cfg['min_samples'] for ss in shared.values()):
                    category = 'ALTERNATIVE_JUNCTION_SUPPORTED_CONNECTION_UNRESOLVED'
                primary = sets['exact']
                primary_label = 'exact'
                r.update(left_reads=len(sets['left']), right_reads=len(sets['right']), coexpressed_samples=len(coexp),
                         bridge_reads=len(sets['bridge']), gene_context_bridge_reads=len(sets['gene_context_bridge']), exact_reads=len(sets['exact']),
                         best_alternative_reads=max([count(s) for s in alt.values()]+[(0,0)])[0],
                         best_alternative_samples=max([count(s) for s in alt.values()]+[(0,0)])[1],
                         best_shared_boundary_reads=max([count(s) for s in shared.values()]+[(0,0)])[0],
                         best_shared_boundary_samples=max([count(s) for s in shared.values()]+[(0,0)])[1])
            else:
                primary = sets['support']
                primary_label = 'support'
                n, ns = count(primary)
                category = 'INSUFFICIENT_SUPPORT'
                if n >= cfg['min_reads'] and ns >= cfg['min_samples']:
                    category = 'CROSS_GENE_CONNECTION_REVIEW' if event['kind'] == 'cross_gene_junction' else 'NOVEL_EXON_CONTEXT_REVIEW'
                    if event['kind'] == 'cross_gene_junction' and count(sets['preserve_both'])[0] >= cfg['min_reads'] and count(sets['preserve_both'])[1] >= cfg['min_samples']:
                        category = 'PRESERVES_BOTH_MODELS_READTHROUGH_COMPATIBLE'
            contexts = event.get('context', [[event['start'], event['end']]])
            motifs = [motif(fa, event['chrom'], a, b, event['strand']) for a, b in contexts]
            canonical = all(m in CANONICAL for m in motifs)
            if not canonical and event['kind'] != 'reference_junction' and category != 'INSUFFICIENT_SUPPORT':
                category = 'NONCANONICAL_CONNECTION_REVIEW'
            n, ns = count(primary)
            r.update(category=category, supporting_reads=n, supporting_samples=ns,
                     raw_supporting_segments=raw[(eid, primary_label, '')],
                     preserved_both_reads=len(sets['preserve_both']), preserved_both_samples=count(sets['preserve_both'])[1],
                     motif=';'.join(motifs), canonical=canonical,
                     context=json.dumps(contexts, separators=(',', ':')),
                     supporting_read_definition=cfg['observation_mode']+' deduplicated within sample',
                     recommended_action='Review evidence; no automatic annotation edit')
            for label, mm in p['alternates'].items():
                matches = [m['transcript_id'] for m in mm if model_matches(event, m)]
                r['alternate_'+label] = ';'.join(matches)
            for scope, detail, ss in [('gene_connection', d, ss) for d, ss in alt.items()]+[('shared_boundary_only', d, ss) for d, ss in shared.items()]:
                a, b = map(int, detail.split(':'))
                alt_rows.append(dict(event_id=eid, evidence_scope=scope, start=a, end=b, supporting_reads=len(ss), supporting_samples=count(ss)[1], motif=motif(fa, event['chrom'], a, b, event['strand'])))
                for sid, nn in sorted(Counter(s for s, obs in ss).items()):
                    sample_rows.append(dict(event_id=eid, evidence=scope, detail=detail, sample_id=sid, supporting_reads=nn,
                                            tissue=meta[sid].get('tissue', ''), stage=meta[sid].get('stage', ''), instar=meta[sid].get('instar', '')))
            for label, ss in sets.items():
                for sid, nn in sorted(Counter(s for s, obs in ss).items()):
                    sample_rows.append(dict(event_id=eid, evidence=label, sample_id=sid, supporting_reads=nn,
                                            tissue=meta[sid].get('tissue', ''), stage=meta[sid].get('stage', ''), instar=meta[sid].get('instar', '')))
            rows.append(r)
    fields = ['event_id','kind','genes','chrom','start','end','strand','category','supporting_reads','supporting_samples','raw_supporting_segments','exact_reads','best_alternative_reads','best_alternative_samples','best_shared_boundary_reads','best_shared_boundary_samples','bridge_reads','gene_context_bridge_reads','left_reads','right_reads','coexpressed_samples','preserved_both_reads','preserved_both_samples','motif','canonical','context','supporting_read_definition','recommended_action'] + ['alternate_'+x for x in p['alternates']]
    write_tsv(out/'candidates.tsv', rows, fields)
    write_tsv(out/'sample_support.tsv', sample_rows, ['event_id','evidence','detail','sample_id','supporting_reads','tissue','stage','instar'])
    write_tsv(out/'alternative_connections.tsv', alt_rows, ['event_id','evidence_scope','start','end','supporting_reads','supporting_samples','motif'])
    plots = out/'plots'
    plots.mkdir(exist_ok=True)
    important = [r for r in rows if r['category'] not in ('EXACT_JUNCTION_SUPPORTED', 'INSUFFICIENT_SUPPORT')]
    important.sort(key=lambda r: (r['category'] != 'UNSUPPORTED_HIGH_INFORMATION', -r['supporting_reads'], r['event_id']))
    shown = important[:args.max_plots]
    cards = []
    for r in shown:
        eid = r['event_id']
        cc = Counter({cid: len(ss) for (ee, cid), ss in chain_support.items() if ee == eid})
        svg = plot_event(all_events[eid], p['models'], p['alternates'], chains, cc)
        (plots/(eid+'.svg')).write_text(svg)
        cards.append('<article><h2>'+html.escape(r['genes']+' | '+r['category'])+'</h2><p>'+html.escape(f"{r['supporting_reads']} supporting reads / {r['supporting_samples']} samples for exact event; motif {r['motif']}. Left/right flank reads: {r.get('left_reads','n/a')}/{r.get('right_reads','n/a')}; immediate-flank bridges: {r.get('bridge_reads','n/a')}; gene-context bridges: {r.get('gene_context_bridge_reads','n/a')}; best alternative junction: {r.get('best_alternative_reads','n/a')} reads / {r.get('best_alternative_samples','n/a')} samples; best shared-boundary alternative: {r.get('best_shared_boundary_reads','n/a')} reads / {r.get('best_shared_boundary_samples','n/a')} samples.")+'</p><a href="plots/'+eid+'.svg"><img src="plots/'+eid+'.svg" alt="Review tracks"></a></article>')
    counts = Counter(r['category'] for r in rows)
    page = '''<!doctype html><html lang="en"><meta charset="utf-8"><title>AnnRE review</title><style>body{font:16px system-ui;max-width:1300px;margin:auto;padding:24px;color:#203040}article{border-top:1px solid #aaa;margin-top:30px}img{width:100%}h2{font-size:20px}table{border-collapse:collapse}td,th{border:1px solid #bbb;padding:6px}</style><h1>AnnRE 0.1 alpha — evidence for annotation review</h1><p>Review candidates, not confirmed annotation errors. Reference joins, observed cross-gene junctions, and novel internal exon splice contexts are distinct event types; these counts are not numbers of affected genes.</p><p>Supporting reads are deduplicated observations within sample; PacBio mode groups movie/ZMW, not original RNA molecules. Top-up libraries share sample IDs. Raw segments remain in the audit. “Parents” is not used. Metadata sample IDs do not establish biological independence.</p><p>High information means both flanks are expressed in multiple samples without a passing connecting read. Low information means insufficient evidence to assess a join. Alternative connection support need not validate the original reference boundary. A shared-boundary alternative can support only one portion of a gene; ALTERNATIVE_JUNCTION_SUPPORTED_CONNECTION_UNRESOLVED does not validate the gene connection. Both-model preservation is readthrough-compatible, not a proven gene merge. Repeated noncanonical junctions require separate review.</p><p><a href="candidates.tsv">All candidate/event evidence</a> · <a href="sample_support.tsv">Sample/tissue support</a> · <a href="alternative_connections.tsv">Alternative coordinates</a> · <a href="manifest.json">Run provenance</a></p>'''
    page += '<p>Alternate input QC (unstranded models excluded): '+html.escape(json.dumps(p['alternate_input_qc']))+'</p>'
    page += '<table><tr><th>Category</th><th>Events</th></tr>'+''.join(f'<tr><td>{html.escape(k)}</td><td>{v}</td></tr>' for k,v in sorted(counts.items()))+'</table>'
    page += f'<p>Showing {len(shown)} of {len(important)} review events. Linear genomic scales; displayed alternate models are capped at six/source. Full splice-chain identity ignores transcript-end variation; displayed ends are examples.</p>' + ''.join(cards)+'</html>'
    (out/'index.html').write_text(page)
    dump(out/'manifest.json', dict(version=__version__, prepared_sha256=phash, config=cfg, sources=p['provenance'], alternate_input_qc=p['alternate_input_qc'], code_sha256=p['code_sha256'], libraries=manifests,
                                  event_counts=dict(counts), selected_gene_count=len(cfg.get('genes', [])),
                                  limitations=['No automatic annotation repair', 'Candidate-focused; no whole-genome performance benchmark',
                                               'No independent error labels or precision/recall estimate', 'Contig dictionary match does not establish sequence identity',
                                               'Alignment strand mode requires oriented transcript reads', 'No CDS/protein adjudication',
                                               'Novel exon discovery requires same-gene annotated immediate flanks; incomplete/unassigned contexts are not enumerated']))
    print(json.dumps(dict(counts), indent=2))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--version', action='version', version=__version__)
    sub = ap.add_subparsers(dest='command', required=True)
    prep = sub.add_parser('prepare', help='Validate inputs and freeze a run')
    prep.add_argument('--config', required=True)
    prep.add_argument('--out', required=True)
    prep.set_defaults(func=prepare)
    ex = sub.add_parser('extract', help='Extract one library, suitable for a SLURM array')
    ex.add_argument('--run', required=True)
    ex.add_argument('--library', type=int, required=True)
    ex.set_defaults(func=extract)
    rep = sub.add_parser('report', help='Require all libraries, deduplicate samples and generate review outputs')
    rep.add_argument('--run', required=True)
    rep.add_argument('--max-plots', type=int, default=30)
    rep.set_defaults(func=report)
    args = ap.parse_args()
    try:
        args.func(args)
    except (ValueError, OSError, KeyError, IndexError) as error:
        ap.exit(2, f'annre: {error}\n')

if __name__ == '__main__':
    main()
