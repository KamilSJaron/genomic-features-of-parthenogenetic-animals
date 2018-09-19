##### GENERATE ANNOTATION FILE UNDERSTANDABLE BY
# gff3
# Avag1
module add UHTS/Assembler/cufflinks/2.2.1;

mkdir -p data/$SP/MCScanX
cp /scratch/beegfs/monthly/kjaron/review-of-asexual-genomes/data/$SP/annotation* data/$SP
cp /scratch/beegfs/monthly/kjaron/review-of-asexual-genomes/data/$SP/genome.fa.gz data/$SP

# scripts/download_data.sh $SP genome
GENOME=data/$SP/genome.fa.gz
# scripts/download_data.sh $SP annotation
GFF=$(echo data/$SP/annotation*)
ll -h $GENOME $GFF
MCScanX_DIR=data/$SP/MCScanX

GENOME_UNZIPED=${GENOME%.*}
GFF_UNZIPED=${GFF%.*}
PROTEINS=$MCScanX_DIR/annotation_proteins.fa

zcat $GENOME > $GENOME_UNZIPED
# sed 's/gi.*.ref.//g' $GENOME_UNZIPED | sed 's/|//g' > data/Pfor1/genome_corrected.fa
# rotifers
# awk '/>/{ print ">"$7 } !/>/ { print $0 } ' $GENOME_UNZIPED | sed 's/,//g' > data/$SP/genome_corrected.fa
# plectus sambesii
# awk '/>/{ print ">PSAMB."$6 } !/>/{print $0} ' $GENOME_UNZIPED | sed s/,$// > data/$SP/genome_corrected.fa
# GENOME_UNZIPED=data/$SP/genome_corrected.fa
zcat $GFF > $GFF_UNZIPED
# sed -i'' -e 's/;Name=;Name=/;Name=/' $GFF_UNZIPED
# Ps591
# sed -i'' -e 's/|/_/' $GFF_UNZIPED
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
}' $GFF_UNZIPED > data/$SP/MCScanX/"$SP"_prot.gff

# Minc1
awk '($3 == "mRNA") {
    OFS="\t";
    match($9, /ID=transcript:Minc3s[0-9]+[a-z]+[0-9]+/);
    print $1, substr($9,RSTART+3,RLENGTH-3), $4, $5
}' $GFF_UNZIPED > data/$SP/MCScanX/"$SP"_prot.gff

# lcla, Mjav1, Mjav2, Mare1, Mare2, Obir1
awk '($3 == "mRNA") {
    OFS="\t";
    match($9, /ID=.+;Parent/);
    print $1, substr($9,RSTART+3,RLENGTH-10), $4, $5
}' $GFF_UNZIPED > $MCScanX_DIR/"$SP"_prot.gff

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

# Minc2 >g19756.t1
awk '($3 == "mRNA") {
    OFS="\t";
    match($9, /ID=[a-z]+[0-9]+.t1/);
    print $1, substr($9,RSTART+3,RLENGTH-3), $4, $5
}' $GFF_UNZIPED > data/$SP/MCScanX/"$SP"_prot.gff

TARGET_DIR=/scratch/beegfs/monthly/kjaron/review-of-asexual-genomes/data/$SP/MCScanX/
mkdir -p $TARGET_DIR
cp data/$SP/MCScanX/"$SP"_prot.gff $TARGET_DIR
cp $PROTEINS $TARGET_DIR