import random
import sys

import gfflib
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import fastalib as fl

MARGIN_SIZE = 150
DONOR = 'GT'
ACCEPTOR = 'AG'

shroom_name = sys.argv[1]
assembly_fasta = sys.argv[2]
gff_file = sys.argv[3]
output_folder = sys.argv[4]

# ---------------------------
# shroom_name = 'Aciaci1'
# gff_file = f'/home/anhvu/Desktop/ascomycota-data/GFFS/{shroom_name}_GeneCatalog_genes_20160228.gff'
# assembly_fasta = f'/home/anhvu/Desktop/mykointrons-data/data/Assembly/{shroom_name}_AssemblyScaffolds.fasta'

db = gfflib.parse_gff(gff_file)

intergenic_pseudo_donors = []
intergenic_pseudo_acceptors = []

with open(assembly_fasta, 'r') as f:
    scaffolds_input = {desc: seq for desc, seq in fl.read_fasta(f)}

for gene in db.all_features(featuretype='gene'):

    if random.randint(0, 10) < 9:
        continue

    for f in db.children(gene['name'][0], featuretype='exon'):
        # Extract the exon described in the GFF file from the given FASTA.
        exon_seq = f.sequence(assembly_fasta)

        scaffold = f.chrom
        strand = f.strand

        cropped_exon = exon_seq[MARGIN_SIZE:-MARGIN_SIZE]

        for position, dimer in fl.dimers(cropped_exon):
            if dimer == DONOR or dimer == ACCEPTOR:
                # compute start and end for the window
                win_start = position - MARGIN_SIZE - 1
                win_end = position + MARGIN_SIZE + 1

                window = cropped_exon[win_start:win_end]

                if win_start < 0 or len(window) != 2 * MARGIN_SIZE + 2:
                    continue    # Insufficient lengths... Skipping the candidate

                assert window in exon_seq

                real_start = f.start - 1 + MARGIN_SIZE + win_start if strand == '+' else f.end - MARGIN_SIZE - win_end
                real_end = f.start - 1 + MARGIN_SIZE + win_end if strand == '+' else f.end - MARGIN_SIZE - win_start

                if strand == '+':
                    assert scaffolds_input[scaffold][real_start:real_end] == window
                if strand == '-':
                    assert Seq(scaffolds_input[scaffold][real_start:real_end]).reverse_complement() == window

                pseudosplicesite = SeqRecord(
                    id=f'{scaffold} {strand} {real_start} {real_end}',
                    seq=Seq(window)
                )

                if dimer == DONOR:
                    intergenic_pseudo_donors.append(pseudosplicesite)
                else:
                    intergenic_pseudo_acceptors.append(pseudosplicesite)

with open(f'{shroom_name}-donor-false.fasta', 'w') as f:
    SeqIO.write(intergenic_pseudo_donors, f, 'fasta')

with open(f'{shroom_name}-acceptor-false.fasta', 'w') as f:
    SeqIO.write(intergenic_pseudo_acceptors, f, 'fasta')
