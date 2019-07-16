### RAW READS

KMER=21

ls raw_reads/*.fastq.gz > FILES
mkdir -p tmp
kmc -k21 -t8 -m32 -ci1 -cs50000 @FILES kmc_kmer_counts tmp

HIST=kmc_$(head -1 FILES  | cut -f 2 -d / | cut -f 1 -d _)_k"$KMER".hist
kmc_tools transform kmc_kmer_counts histogram $HIST -cx50000

# kmc -k17 -t60 -m64 -ci1 -cs500000000 @TRIMMED_FILES pvir1_kmc_trimmed_kmer_counts tmp
# kmc_tools transform pvir1_kmc_trimmed_kmer_counts histogram temp.hist -cx500000000
# grep -v "[[:space:]]0" temp.hist > kmc_trimmed_$HIST_SUFFIX # cleaning up the rows with 0 coverage

ls trimmed_reads/*.fastq.gz > FILES
mkdir -p tmp logs kmc
bsub -n 8 -M 32000000 -o logs/kmc_counts_trimmed.log "kmc -k21 -t8 -m32 -ci1 -cs50000 @FILES kmc/kmc_kmer_counts tmp"

bsub -o logs/kmc_hist.log "kmc_tools transform kmc/kmc_kmer_counts histogram temp.hist -cx50000"

KMER=21
HIST=kmc/kmc_$(head -1 FILES  | cut -f 2 -d / | cut -f 1 -d -)_k"$KMER".hist
grep -v "[[:space:]]0" temp.hist > $HIST