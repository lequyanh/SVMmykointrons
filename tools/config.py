ASSEMBLIES_LOC = "/home/anhvu/Desktop/mykointrons-data/data/Assembly"
NEWSEQUENCES_LOC = "/home/anhvu/Desktop/mykointrons-data/new-sequences"
GFFS_LOC = '/home/anhvu/Desktop/mykointrons-data/data/GFF'

EXON_CSV_SUFF = "_exon_positions.csv"
EXON_FASTA_SUFF = "_exons.fasta"

INTRON_CSV_SUFF = "_intron_positions.csv"
INTRON_FASTA_SUFF = "-introns.fasta"


def get_fungi_assembly(fungi_name: str):
    assembly_fasta = f'{ASSEMBLIES_LOC}/{fungi_name}_AssemblyScaffolds.fasta'
    return assembly_fasta


def get_fungi_intron_fasta(fungi_name: str):
    introns_fasta = f'{NEWSEQUENCES_LOC}/{fungi_name}/{fungi_name}-introns.fasta'
    # Uncomment if reduced assemblies and introns is to use
    # introns_fasta = f'/home/anhvu/Desktop/mykointrons-data/reduced/new-sequences/{fungi_name}-introns.fasta'
    return introns_fasta


def get_fungi_exons_fasta(fungi_name: str):
    exons_file = f'{NEWSEQUENCES_LOC}/{fungi_name}/{fungi_name}_exons.fasta'
    return exons_file


def get_fungi_exons_positions(fungi_name: str):
    exons_positions = f'{NEWSEQUENCES_LOC}/{fungi_name}/{fungi_name}_exon_positions.csv'
    return exons_positions


def get_fungi_false_wins_fasta(fungi_name: str, strand: str, site: str):
    suffix = 'false-intragenic'
    return _get_fungi_wins_fasta(fungi_name, strand, site, suffix)


def get_fungi_true_wins_fasta(fungi_name: str, strand: str, site: str):
    suffix = 'true'
    return _get_fungi_wins_fasta(fungi_name, strand, site, suffix)


def _get_fungi_wins_fasta(fungi_name: str, strand: str, site: str, suffix: str):
    strand_suff = 'plus' if strand == '+' else 'minus'
    return f'{NEWSEQUENCES_LOC}/{fungi_name}/{fungi_name}-{site}-{suffix}-{strand_suff}.fasta'
