#!/bin/bash
# 1. argument should be name of stick insect (5_Tge for instance)
# 2. argument is insert size

SP="$1"
REF="$2"

BAM=map_to_"$REF".bam
PATH_TO_DATA=/scratch/beegfs/monthly/kjaron/review-of-asexual-genomes/data/"$SP"/


if [[ -s "$PATH_TO_DATA"/"$READS1" ]]
then
	echo "$PATH_TO_DATA"/"$READS1" ALREADY EXISTS
	exit 1
fi

if [[ ! -s "$PATH_TO_DATA"/"$BAM" ]]
then
	echo "$PATH_TO_DATA"/"$BAM" IS MISSING
	exit 1
fi

bsub <<< """
#BSUB -L /bin/bash
#BSUB -J "$SP"resc
#BSUB -q bgee
#BSUB -o "$SP"_rescue.out
#BSUB -e "$SP"_rescue.err
#BSUB -n 16
#BSUB -M 24000000
#BSUB -R \"rusage[tmp=120000]\"

module add UHTS/Analysis/BEDTools/2.26.0;
module add UHTS/Analysis/samtools/1.3;

# make local directory for computations
LOCAL_DIR=/scratch/local/monthly/kjaron/"$SP"_"$REF"_rescue
mkdir -p \$LOCAL_DIR/temp
export TMPDIR=\$LOCAL_DIR/temp

cd \$LOCAL_DIR
cp "$PATH_TO_DATA"/"$BAM" .

# filter duplicates and sort bam by name
samtools view -h -F 256 "$BAM" | samtools sort -@15 -n - -o aln.sorted.bam

# convert sorted bam with unique reads to fq
bedtools bamtofastq -i aln.sorted.bam \
                      -fq reads_R1.fq \
                      -fq2 reads_R2.fq

# zip reads
pigz -p 16 reads_R1.fq
pigz -p 16 reads_R2.fq

# move zipped reads back
mkdir -p "$PATH_TO_DATA"

# move to working dir
mv reads_R1.fq.gz "$PATH_TO_DATA" &
mv reads_R2.fq.gz "$PATH_TO_DATA" &

wait

rm "$BAM" aln.sorted.bam
rm -r temp
rmdir \$LOCAL_DIR

"""





