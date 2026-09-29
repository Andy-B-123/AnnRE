"""Coordinate-safe, annotation-independent evidence primitives (0-based half-open)."""
import hashlib
import re
from collections import defaultdict


def ident(value):
    return hashlib.sha256(str(value).encode()).hexdigest()[:20]


def attributes(text):
    # Semicolons inside quoted GTF values are not separators.
    return dict(re.findall(r'(\S+)\s+"((?:[^"\\]|\\.)*)"\s*;?', text))


def load_gtf(path, allow_unstranded=False, stats=None):
    models = {}
    unstranded = set()
    with open(path) as handle:
        for line in handle:
            if line.startswith('#') or not line.strip():
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) != 9:
                raise ValueError(f'Invalid GTF line: {line[:80]}')
            if f[2] != 'exon':
                continue
            a = attributes(f[8])
            tid, gid = a.get('transcript_id'), a.get('gene_id')
            if not tid or not gid:
                raise ValueError('Exons require gene_id and transcript_id')
            if f[6] == '.' and allow_unstranded:
                unstranded.add(tid)
                continue
            if f[6] not in ('+', '-'):
                raise ValueError('Exons require +/- strand')
            start, end = int(f[3])-1, int(f[4])
            if start < 0 or end <= start:
                raise ValueError('Invalid exon coordinates')
            m = models.setdefault(tid, dict(transcript_id=tid, gene_id=gid,
                                            chrom=f[0], strand=f[6], exons=[]))
            if (m['gene_id'], m['chrom'], m['strand']) != (gid, f[0], f[6]):
                raise ValueError(f'Inconsistent transcript: {tid}')
            m['exons'].append((start, end))
    if stats is not None:
        stats['excluded_unstranded_transcripts'] = len(unstranded)
        stats['included_transcripts'] = len(models)
    if not models and not (allow_unstranded and unstranded):
        raise ValueError('No exon models in GTF (GFF3/CDS-only input is not supported)')
    for m in models.values():
        m['exons'] = sorted(set(m['exons']))
        if any(a[1] >= b[0] for a, b in zip(m['exons'], m['exons'][1:])):
            raise ValueError(f'Overlapping/adjacent exons in {m["transcript_id"]}')
        m['chain'] = [(a[1], b[0]) for a, b in zip(m['exons'], m['exons'][1:])]
    return models


def geometry(start, cigar):
    blocks, gaps, exons = [], [], []
    pos = exon_start = start
    for op, length in cigar or []:
        if op in (0, 7, 8):
            blocks.append((pos, pos+length))
            pos += length
        elif op == 2:
            pos += length
        elif op == 3:
            exons.append((exon_start, pos))
            gaps.append((pos, pos+length))
            pos += length
            exon_start = pos
    exons.append((exon_start, pos))
    return blocks, gaps, exons


def overlap(blocks, interval):
    s, e = interval
    return sum(max(0, min(b, e)-max(a, s)) for a, b in blocks)


def union(intervals):
    result = []
    for s, e in sorted(intervals):
        if result and s <= result[-1][1]:
            result[-1] = (result[-1][0], max(e, result[-1][1]))
        else:
            result.append((s, e))
    return result


def observation(name, mode):
    parts = name.split('/')
    if mode == 'pacbio-zmw':
        if len(parts) < 3 or not parts[1].isdigit():
            raise ValueError('Unrecognised PacBio movie/ZMW read identifier')
        return '/'.join(parts[:2])
    return name


def rejection(read, config):
    flags = []
    for attr in ('is_unmapped', 'is_secondary', 'is_supplementary', 'is_duplicate', 'is_qcfail'):
        if getattr(read, attr):
            flags.append(attr)
    if read.has_tag('SA'):
        flags.append('SA')
    if read.mapping_quality == 255 or read.mapping_quality < config['min_mapq']:
        flags.append('mapq')
    clip = sum(n for op, n in read.cigartuples or [] if op in (4, 5))
    if clip / max(1, (read.query_length or 0) + sum(n for op, n in read.cigartuples or [] if op == 5)) > config['max_clip_fraction']:
        flags.append('clipping')
    if not read.has_tag('NM') or read.get_tag('NM') / max(1, read.query_alignment_length or 0) > config['max_edit_fraction']:
        flags.append('NM')
    return flags


def motif(fasta, chrom, start, end, strand):
    left, right = fasta.fetch(chrom, start, start+2).upper(), fasta.fetch(chrom, end-2, end).upper()
    if strand == '-':
        rc = lambda s: s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
        left, right = rc(right), rc(left)
    return left+'-'+right


CANONICAL = {'GT-AG', 'GC-AG', 'AT-AC'}


def classify_reference(exact, alternatives, bridges, left, right, coexpressed, cfg):
    enough = lambda x: x[0] >= cfg['min_reads'] and x[1] >= cfg['min_samples']
    if enough(exact):
        return 'EXACT_JUNCTION_SUPPORTED'
    if any(enough(x) for x in alternatives):
        return 'ALTERNATIVE_CONNECTION_SUPPORTED'
    if exact[0] or bridges[0]:
        return 'CONNECTION_EVIDENCE_REVIEW'
    if min(left, right) >= cfg['min_flank_reads'] and coexpressed >= cfg['min_samples']:
        return 'UNSUPPORTED_HIGH_INFORMATION'
    return 'UNSUPPORTED_LOW_INFORMATION'


def subchain(short, long):
    needle = [tuple(x) for x in short]
    haystack = [tuple(x) for x in long]
    return bool(needle) and any(needle == haystack[i:i+len(needle)] for i in range(len(haystack)-len(needle)+1))


def preserved(model, gaps, blocks, anchor):
    ex = model['exons']
    if not model['chain']:
        return overlap(blocks, ex[0]) >= .95*(ex[0][1]-ex[0][0])
    return subchain(model['chain'], gaps) and all(overlap(blocks, e) >= min(anchor, e[1]-e[0]) for e in (ex[0], ex[-1]))
