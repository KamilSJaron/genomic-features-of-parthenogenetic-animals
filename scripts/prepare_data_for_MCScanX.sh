##### GENERATE ANNOTATION FILE UNDERSTANDABLE BY
# gff3
# Avag1
GENOME_UNZIPED=${GENOME%.*}
GFF_UNZIPED=${GFF%.*}
PROTEINS=data/$SP/$(basename $GFF .gff3.gz)_proteins.fa

zcat $GENOME > $GENOME_UNZIPED
zcat $GFF > $GFF_UNZIPED

gffread -g $GENOME_UNZIPED -y $PROTEINS $GFF_UNZIPED
sed -i'' -e 's/.$//' $PROTEINS

awk '($3 == "mRNA") {
    OFS="\t";
    match($9, /Name=.+/);
    print $1, substr($9,RSTART+5,RLENGTH), $4, $5
}' "$GFF_UNZIPED" > "$MCScanX_DIR"/"$SP"_prot.gff

# gff - Dpac
awk '($3 == "mRNA") {
    OFS="\t";
    match($9, /ID=[a-z]+[0-9]+/);
    print $1, substr($9,RSTART+3,RLENGTH-3), $4, $5
}' annotation.gff > MCScanX/Dpac1_prot.gff
MCScanX MCScanX/Dpac1_prot

# lcla
awk '($3 == "mRNA") {
    OFS="\t";
    match($9, /ID=.+;Parent/);
    print $1, substr($9,RSTART+3,RLENGTH-10), $4, $5
}' lclav_annotation.gff3.1 > lclav_prot.gff