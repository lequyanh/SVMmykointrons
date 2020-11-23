############################
# FUNCTION FOR DOWNLOADING #
############################

function gdrive_download () {
 CONFIRM=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate "https://docs.google.com/uc?export=download&id=$1" -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')
 wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$CONFIRM&id=$1" -O $2
 rm -rf /tmp/cookies.txt
}

#################
# PRUNE CONTIGS #
#################

# input files
metagenome_file=$1
cut_fwd_file=$2
cut_rev_file=$3

# get affected contigs
awk -F';' '{print ">"$1" "}' $cut_fwd_file > fwd_ctgs_names.txt
awk -F';' '{print ">"$1" "}' $cut_rev_file > rev_ctgs_names.txt

grep --no-group-separator -A 1 -F -f fwd_ctgs_names.txt $metagenome_file > ${metagenome_file}.fwd
grep --no-group-separator -A 1 -F -f rev_ctgs_names.txt $metagenome_file > ${metagenome_file}.rev

#gdrive_download 1_-gzSgegMvIinpAZHNS2iJsSn1wzqeK5 prune_contigs.py

python2.7 prune_contigs.py ${metagenome_file}.fwd $cut_fwd_file  ${metagenome_file}.fwd.cut
python2.7 prune_contigs.py ${metagenome_file}.rev $cut_rev_file  ${metagenome_file}.rev.cut

perl -p -e "s/^>/>ori-/g" $metagenome_file > ori_to_mg_rast.fa
perl -p -e "s/^>/>fwd-/g" ${metagenome_file}.fwd.cut > fwd_to_mg_rast.fa
perl -p -e "s/^>/>rev-/g" ${metagenome_file}.rev.cut > rev_to_mg_rast.fa

cat ori_to_mg_rast.fa fwd_to_mg_rast.fa rev_to_mg_rast.fa > all_to_mg_rast.fa

rm ori_to_mg_rast.fa
rm fwd_to_mg_rast.fa
rm rev_to_mg_rast.fa

# all_to_mg_rast.fa - will be annotated by (mg-rast or "CUSTOM ANNOTATION")...

###########
# MG-RAST #
###########

# DOWNLOAD
protein_file=mgm4908545.3.350.genecalling.coding.faa
tax_file=taxonomy_mgm4908545.3_img.csv

#export LANG=en_US.UTF-8
#export LC_ALL=en_US.UTF-8
#sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr $tax_file | sort -u -k1,1 --merge > ${tax_file}_best.txt
#
## format annotation
## mgm4912149.3|fwd_k141_10000242_3854_6888_+|IMG  8a1b5566beaf0f33989c29ff576d4629        56.96   79      34              879     957     388     466     1.3e-19 94.0    [Mycobacterium kansasii ATCC 12478]
## format blast output - out 6
## 1- qseqid	 2.- sseqid	 3.- pident	 4.- length	 5.- mismatch	 6.- gapopen	 7.- qstart	 8.- qend	 9.- sstart	 10.- send	 11.- evalue	 12.- bitscore
#awk -F'[|\t]' '{print $2 "\t" $15 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14}' ${tax_file}_best.txt > ${tax_file}_best_reformate.txt

tax_file_best=${tax_file}_best_reformate.txt

########################
## "CUSTOM ANNOTATION" #
########################
#
# #custom gene call or use of MG-RAST genecall ($protein_file)
# perl /home/ubuntu/FGS-master/run_FragGeneScan.pl -complete 0 -train illumina_10 -genome all_to_mg_rast.fa -out genecalling_custom.faa
#
# protein_file=genecalling_custom.faa
#
# # BLAST JGI FUNGAL CUSTOM DATABASE
# diamond blastp -d /mnt/DATA/JGI_FUNGAL_PROTEINS_ANNOTATED_20200511 -q $protein_file -e 1E-5 -o ${protein_file}.out6.txt -f 6 -p 20
# sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr ${protein_file}.out6.txt | sort -u -k1,1 --merge > ${protein_file}.out6_best.txt
#
# # formate annotation
# gdrive_download 1TKLRzGehrNfaXteURVKPt_3-7e5rleCw jgi_abr_org_list.txt
# gdrive_download 1XBTtiC1JYl2rzeV7idN2WrveEZknmnQi replace_fungal_annot_by_taxname.py
# python2.7 replace_fungal_annot_by_taxname.py ${protein_file}.out6_best.txt jgi_abr_org_list.txt ${protein_file}.out6_best_reformate.txt
#
##################################
## compare different annotations #
##################################
#
#gdrive_download 11Iwd16M7RUUC56CiNTcTCbUskJVB37i8 get_best_annotation_by_biscore.py
python2.7 get_best_annotation_by_biscore.py taxonomy_mgm4908545.3_img.csv_best_reformate.txt taxonomy_mgm4908545.3_fungi_best_reformate.txt best_of_annotations_best.txt

tax_file_best=best_of_annotations_best.txt
#
####################
# FINAL COMPARISON #
####################

grep 'ori-' $tax_file_best > ori_final_annot.txt
grep 'fwd-' $tax_file_best > fwd_final_annot.txt
grep 'rev-' $tax_file_best > rev_final_annot.txt

#gdrive_download 1wLTAnEFr9_21XLkMifULNS7XmaOEuZMy compare_annotation_final.py

python2.7 compare_annotation_final.py ori_final_annot.txt fwd_final_annot.txt rev_final_annot.txt compared_annotatoin.txt


###############
# GET RESULTS #
###############

#awk -F'\t' '{print $12}' compared_annotatoin.txt | sort | uniq -c > compared_annotatoin_random_bitscore_sums.txt
#awk -F'\t' '{print $12}' ori_final_annot.txt | sort | uniq -c > original_annotatoin_random_bitscore_sums.txt

