import argparse
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path
import pysam
from annre.core import (attributes, geometry, observation, motif, rejection, classify_reference, load_gtf)
from annre.cli import prepare, extract, report, read_tsv


class Primitives(unittest.TestCase):
    def test_gtf_quoted_semicolon(self):
        self.assertEqual(attributes('gene_id "g"; note "one; two"; transcript_id "t";')['note'], 'one; two')

    def test_unstranded_alternate_explicit(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'alt.gtf'
            path.write_text('chr\ttest\texon\t1\t100\t.\t.\t.\tgene_id "g"; transcript_id "t";\n')
            with self.assertRaisesRegex(ValueError, 'strand'):
                load_gtf(path)
            stats = {}
            self.assertEqual(load_gtf(path, allow_unstranded=True, stats=stats), {})
            self.assertEqual(stats['excluded_unstranded_transcripts'], 1)

    def test_deletion_not_intron(self):
        blocks, gaps, exons = geometry(100, [(0,30),(2,10),(0,20),(3,100),(0,40)])
        self.assertEqual(gaps, [(160,260)])
        self.assertEqual(exons, [(100,160),(260,300)])
        self.assertEqual(blocks, [(100,130),(140,160),(260,300)])

    def test_zmw(self):
        self.assertEqual(observation('movie/42/1', 'pacbio-zmw'), observation('movie/42/2', 'pacbio-zmw'))
        with self.assertRaises(ValueError):
            observation('read42', 'pacbio-zmw')
        self.assertEqual(observation('ont-id', 'read-name'), 'ont-id')

    def test_quality(self):
        read = pysam.AlignedSegment()
        read.query_name = 'm/1/1'
        read.query_sequence = 'A'*100
        read.cigartuples = [(0,100)]
        read.mapping_quality = 255
        cfg = dict(min_mapq=20,max_clip_fraction=.1,max_edit_fraction=.05)
        self.assertIn('mapq', rejection(read,cfg))
        self.assertIn('NM', rejection(read,cfg))
        read.mapping_quality = 60
        read.set_tag('NM', 0)
        self.assertEqual(rejection(read,cfg), [])
        read.set_tag('SA','chr,1,+,100M,60,0;')
        self.assertIn('SA', rejection(read,cfg))

    def test_alternative_context_not_pooled(self):
        cfg = dict(min_reads=3,min_samples=2,min_flank_reads=10)
        self.assertEqual(classify_reference((0,0),[(2,2),(2,2)],(4,2),30,40,2,cfg), 'CONNECTION_EVIDENCE_REVIEW')


class Workflow(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.fa = self.root/'genome.fa'
        seq = list('A'*7000)
        for a,b in [(200,400),(1100,1200),(2100,2200),(2300,2400),(2500,2600),(4100,4200),(4300,4400),(5100,5200)]:
            seq[a:a+2] = 'GT'
            seq[b-2:b] = 'AG'
        # Negative-strand GT-AG intron: genomic CT ... AC.
        seq[6000:6002], seq[6098:6100] = 'CT', 'AC'
        self.fa.write_text('>chr\n'+''.join(seq)+'\n')
        pysam.faidx(str(self.fa))
        models = {'split':[(100,200),(400,500)], 'supported':[(1000,1100),(1200,1300)],
                  'left':[(2000,2100),(2200,2300)], 'right':[(2400,2500),(2600,2700)],
                  'novel':[(4000,4100),(4400,4500)], 'low':[(5000,5100),(5200,5300)],
                  'alternative':[(5500,5600),(5700,5800),(6200,6300)],
                  'shared':[(6500,6600),(6800,6900)]}
        gtf = self.root/'ref.gtf'
        gtf.write_text(''.join(f'chr\ttest\texon\t{a+1}\t{b}\t.\t+\t.\tgene_id "{g}"; transcript_id "{g}.1";\n' for g,ex in models.items() for a,b in ex))
        entries = []
        # Topup sample A deliberately repeats ZMW IDs from its first library.
        for lib, sample in [('a','A'),('b','B'),('topup','A')]:
            bam = self.root/(lib+'.bam')
            with pysam.AlignmentFile(str(bam),'wb',header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chr','LN':7000}]}) as h:
                for start,cigar,n,base in [(100,[(0,100)],12,1),(400,[(0,100)],12,100),
                       (1000,[(0,100),(3,100),(0,100)],4,200),
                       (2000,[(0,100),(3,100),(0,100),(3,100),(0,100),(3,100),(0,100)],4,300),
                       (4000,[(0,100),(3,100),(0,100),(3,100),(0,100)],4,400),
                       (5500,[(0,100),(3,200),(0,100),(3,300),(0,100)],4,500),
                       (6670,[(0,30),(3,100),(0,100)],4,600)]:
                    for k in range(n):
                        r = pysam.AlignedSegment()
                        r.query_name = f'{sample}/{base+k}/1'
                        r.query_sequence = 'A'*sum(l for op,l in cigar if op == 0)
                        r.flag = 0
                        r.reference_id = 0
                        r.reference_start = start
                        r.mapping_quality = 60
                        r.cigartuples = cigar
                        r.set_tag('NM',0)
                        h.write(r)
            pysam.index(str(bam))
            entries.append(f'{lib}\t{sample}\t{bam}\tTissue\tStage\n')
        samples = self.root/'samples.tsv'
        samples.write_text('library_id\tsample_id\tbam\ttissue\tstage\n'+''.join(entries))
        config = self.root/'config.json'
        config.write_text(json.dumps(dict(genome=str(self.fa),reference=str(gtf),samples=str(samples),context_bp=0,alternates={'same':str(gtf)})))
        self.run = self.root/'run'
        with contextlib.redirect_stdout(io.StringIO()):
            prepare(argparse.Namespace(config=str(config),out=str(self.run)))

    def tearDown(self):
        self.tmp.cleanup()

    def test_minus_motif(self):
        with pysam.FastaFile(str(self.fa)) as fa:
            self.assertEqual(motif(fa,'chr',6000,6100,'-'), 'GT-AG')

    def test_incomplete_array_refused(self):
        with self.assertRaisesRegex(ValueError,'Incomplete run'):
            report(argparse.Namespace(run=str(self.run),max_plots=5))

    def test_end_to_end(self):
        with contextlib.redirect_stdout(io.StringIO()):
            for i in range(3):
                extract(argparse.Namespace(run=str(self.run),library=i))
            report(argparse.Namespace(run=str(self.run),max_plots=10))
        rows = read_tsv(self.run/'candidates.tsv')
        ref = {r['genes']:r for r in rows if r['kind'] == 'reference_junction'}
        self.assertEqual(ref['split']['category'],'UNSUPPORTED_HIGH_INFORMATION')
        self.assertEqual(ref['low']['category'],'UNSUPPORTED_LOW_INFORMATION')
        shifted = next(r for r in rows if r['kind'] == 'reference_junction' and r['start'] == '5800')
        self.assertEqual(shifted['category'], 'ALTERNATIVE_CONNECTION_SUPPORTED')
        self.assertEqual(shifted['bridge_reads'], '0')
        self.assertEqual(shifted['gene_context_bridge_reads'], '8')
        self.assertEqual(ref['shared']['category'],'ALTERNATIVE_JUNCTION_SUPPORTED_CONNECTION_UNRESOLVED')
        self.assertEqual(ref['shared']['best_shared_boundary_reads'],'8')
        self.assertEqual(ref['shared']['gene_context_bridge_reads'],'0')
        self.assertEqual(ref['supported']['category'],'EXACT_JUNCTION_SUPPORTED')
        self.assertEqual(ref['supported']['supporting_reads'],'8')
        self.assertEqual(ref['supported']['raw_supporting_segments'],'12')
        self.assertEqual(ref['supported']['supporting_samples'],'2')
        cross = [r for r in rows if r['kind']=='cross_gene_junction']
        self.assertEqual(len(cross),1)
        self.assertEqual(cross[0]['category'],'PRESERVES_BOTH_MODELS_READTHROUGH_COMPATIBLE')
        novel = [r for r in rows if r['kind']=='novel_exon_context' and r['genes']=='novel']
        self.assertEqual(len(novel),1)
        self.assertEqual(novel[0]['category'],'NOVEL_EXON_CONTEXT_REVIEW')
        self.assertEqual(novel[0]['supporting_reads'],'8')
        self.assertEqual(ref['supported']['alternate_same'],'supported.1')
        self.assertTrue((self.run/'index.html').exists())
        # Tampering and incomplete evidence must not be silently consumed.
        (self.run/'evidence/0000/events.json').write_text('{}')
        with self.assertRaisesRegex(ValueError,'checksum mismatch'):
            report(argparse.Namespace(run=str(self.run),max_plots=5))

if __name__ == '__main__':
    unittest.main()
