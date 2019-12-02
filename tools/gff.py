import random

import gffutils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import fastalib as fl

MARGIN_SIZE = 150
DONOR = 'GT'
ACCEPTOR = 'AG'

fn = gffutils.example_filename('/home/anhvu/Desktop/ascomycota-data/GFFS/Aciaci1_GeneCatalog_genes_20160228.gff')
assembly_fasta = '/home/anhvu/Desktop/mykointrons-data/data/Assembly/Aciaci1_AssemblyScaffolds.fasta'


def transform_func(d):
    try:
        d['transcriptId'] = d['proteinId']
    except KeyError:
        pass
    return d


db = gffutils.create_db(
    fn, ":memory:",
    id_spec={'transcript': 'transcriptId', 'gene': 'name'},
    gtf_transcript_key='transcriptId',
    gtf_gene_key='name',
    transform=transform_func
)

intergenic_pseudo_donors = []
intergenic_pseudo_acceptors = []

for gene in db.all_features(featuretype='gene'):

    if random.randint(0, 10) <= 9:
        continue

    for f in db.children(gene['name'][0], featuretype='exon'):
        exon_seq = f.sequence(assembly_fasta)

        scaffold = f.chrom
        strand = f.strand

        cropped_exon = exon_seq[MARGIN_SIZE:-MARGIN_SIZE]

        for position, dimer in fl.dimers(cropped_exon):
            if dimer == DONOR or dimer == ACCEPTOR:
                # compute start and end for the window
                win_start = position - MARGIN_SIZE - 1
                win_end = position + MARGIN_SIZE + 1

                if win_start < 0:
                    print('Donor within margin. Skipping')
                    continue

                window = cropped_exon[win_start:win_end]
                assert window in exon_seq

                real_start = f.start + win_start if strand == '+' else f.end - win_end
                real_end = f.start + win_end if strand == '+' else f.end - win_start

                print(f'{real_start}   -    {real_end}')

                pseudosplicesite = SeqRecord(
                    id=f'{scaffold} {strand} {real_start} {real_end}',
                    seq=Seq(window)
                )

                if dimer == DONOR:
                    intergenic_pseudo_donors.append(pseudosplicesite)
                else:
                    intergenic_pseudo_acceptors.append(pseudosplicesite)

with open('sadfasd-false-donor', 'w') as f:
    SeqIO.write(intergenic_pseudo_donors, f, 'fasta')

with open('sadfasd-false-acceptor', 'w') as f:
    SeqIO.write(intergenic_pseudo_acceptors, f, 'fasta')
