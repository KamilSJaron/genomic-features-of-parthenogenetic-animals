##### GENERATE ANNOTATION FILE UNDERSTANDABLE BY
# gff3
# Avag1
module add UHTS/Assembler/cufflinks/2.2.1;

################
#### README ####
################
# This script contains bits of code that were manually executed for speices mentioned in commnets
# the reason is that evey single annotatoin file is different
# and I had to make sure that I will manage to extract protein sequences that have compatible headars with annotation_proteins
# the two exceptions are Fcan1 and Lcla1 for which I downloaded directly the protein files
################

mkdir -p data/$SP/MCScanX

# scripts/download_data.sh $SP genome
GENOME=data/$SP/genome.fa.gz
# scripts/download_data.sh $SP annotation
GFF=$(echo data/$SP/annotation*)
ll -h $GENOME $GFF
MCScanX_DIR=data/$SP/MCScanX

GENOME_UNZIPED=${GENOME%.*}
GFF_UNZIPED=${GFF%.*}
PROTEINS=data/$SP/annotation_proteins.fa

zcat $GENOME > $GENOME_UNZIPED
# sed 's/gi.*.ref.//g' $GENOME_UNZIPED | sed 's/|//g' > data/Pfor1/genome_corrected.fa
# rotifers A. nanus
# awk '/>/{ print ">"$7 } !/>/ { print $0 } ' $GENOME_UNZIPED | sed 's/,//g' > data/$SP/genome_corrected.fa
# plectus sambesii
# awk '/>/{ print ">PSAMB."$6 } !/>/{print $0} ' $GENOME_UNZIPED | sed s/,$// > data/$SP/genome_corrected.fa
# Dcor1
# awk '/>/{ print ">"$6 } !/>/ { print $0 } ' $GENOME_UNZIPED | sed 's/,//g' > data/$SP/genome_corrected.fa
# GENOME_UNZIPED=data/$SP/genome_corrected.fa
zcat $GFF > $GFF_UNZIPED
# sed -i'' -e 's/;Name=;Name=/;Name=/' $GFF_UNZIPED
# Ps591, Dcor1 as well, A. nanus
# sed -i'' -e 's/|/_/' $GFF_UNZIPED
# sed -i'' -e 's/|size[0-9]*//' $GFF_UNZIPED # Anan1
# and this SEQ={what gffread complains about}
# grep -v $SEQ $GFF_UNZIPED > temp
# rm $GFF_UNZIPED
# mv temp $GFF_UNZIPED
# gffread -g $GENOME_UNZIPED -y $PROTEINS $GFF_UNZIPED

# for P sambesii I had to remove two genes from the annotation (mapping outside of the assembled genome...)

gffread -g $GENOME_UNZIPED -y $PROTEINS $GFF_UNZIPED
sed -i'' -e 's/.$//' $PROTEINS

head -1 $GENOME_UNZIPED
grep ">" $PROTEINS | head
grep "mRNA" $GFF_UNZIPED | head -1

# DEFAULT
# used for lcla, Mjav1, Mjav2, Mare1, Mare2, Obir1
awk '($3 == "mRNA") {
    OFS="\t";
    match($9, /ID=.+;Parent/);
    print $1, substr($9,RSTART+3,RLENGTH-10), $4, $5
}' $GFF_UNZIPED > $MCScanX_DIR/"$SP"_prot.gff

# gff - Dpac, Dcor
awk '($3 == "mRNA") {
    OFS="\t";
    match($9, /ID=.+/);
    print $1, substr($9,RSTART+3,RLENGTH-3), $4, $5
}' $GFF_UNZIPED > data/$SP/MCScanX/"$SP"_prot.gff

# Minc1
awk '($3 == "mRNA") {
    OFS="\t";
    match($9, /ID=transcript:Minc3s[0-9]+[a-z]+[0-9]+/);
    print $1, substr($9,RSTART+3,RLENGTH-3), $4, $5
}' $GFF_UNZIPED > data/$SP/MCScanX/"$SP"_prot.gff

# Tpre1
awk '($3 == "mRNA") {
    OFS="\t";
    match($9, /ID=.+;Name/);
    print $1, substr($9,RSTART+3,RLENGTH-8), $4, $5
}' $GFF_UNZIPED > $MCScanX_DIR/"$SP"_prot.gff

# gff - Pdav
awk '($3 == "transcript") {
    OFS="\t";
    match($9, /ID=[a-z]+[0-9]+.t1/);
    print $1, substr($9,RSTART+3,RLENGTH-3), $4, $5
}' $GFF_UNZIPED > data/$SP/MCScanX/"$SP"_prot.gff

# Minc2
awk '($3 == "mRNA") {
    OFS="\t";
    match($9, /ID=[a-z]+[0-9]+.t1/);
    print $1, substr($9,RSTART+3,RLENGTH-3), $4, $5
}' $GFF_UNZIPED > data/$SP/MCScanX/"$SP"_prot.gff

# Anan
awk '($3 == "transcript") {
    OFS="\t";
    print $1, $9, $4, $5
}' "$GFF_UNZIPED" > "$MCScanX_DIR"/"$SP"_prot.gff

TARGET_DIR=data/$SP/MCScanX/
mkdir -p $TARGET_DIR
cp data/$SP/MCScanX/"$SP"_prot.gff $TARGET_DIR
cp $PROTEINS $TARGET_DIR