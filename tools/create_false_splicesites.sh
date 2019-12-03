SPECIES_LIST="/home/anhvu/PycharmProjects/mycointrons/ascomycotas.txt"
GFF_LOC="/home/anhvu/Desktop/ascomycota-data/GFFS/"
ASSEMBLY_LOC="/home/anhvu/Desktop/mykointrons-data/data/Assembly/"

OUTPUT_LOC="/home/anhvu/Desktop/mykointrons-data/new-sequences/"

while read shroom; do
  echo "${shroom}"
  shroom_assembly="$ASSEMBLY_LOC/${shroom}_AssemblyScaffolds.fasta"

  find $GFF_LOC -name "${shroom}*.gff" -print0 | xargs -I{} python gff.py "${shroom}" "${shroom_assembly}" {} "$OUTPUT_LOC/${shroom}"

  mv "${shroom}-acceptor-false.fasta" $OUTPUT_LOC/"${shroom}"/
  mv "${shroom}-donor-false.fasta" $OUTPUT_LOC/"${shroom}"/

done < $SPECIES_LIST