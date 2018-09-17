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
echo $GENOME $GFF
MCScanX_DIR=data/$SP/MCScanX

GENOME_UNZIPED=${GENOME%.*}
GFF_UNZIPED=${GFF%.*}
PROTEINS=$MCScanX_DIR/annotation_proteins.fa

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
}' $GFF_UNZIPED > data/$SP/MCScanX/"$SP"_prot.gff

# Minc1
awk '($3 == "mRNA") {
    OFS="\t";
    match($9, /ID=transcript:Minc3s[0-9]+[a-z]+[0-9]+/);
    print $1, substr($9,RSTART+3,RLENGTH-3), $4, $5
}' $GFF_UNZIPED > data/$SP/MCScanX/"$SP"_prot.gff

# lcla
awk '($3 == "mRNA") {
    OFS="\t";
    match($9, /ID=.+;Parent/);
    print $1, substr($9,RSTART+3,RLENGTH-10), $4, $5
}' lclav_annotation.gff3.1 > lclav_prot.gff

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