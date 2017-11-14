# from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
# NCBI = NCBIRemoteProvider(email="kamiljaron+ncbi@gmail.com") # email required by NCBI to prevent abuse
# import numpy

species='Cbir Avag Fcan Lcla Dcor Dpac Minc1 Minc2 Mjav Mare Mflo Mhap'
# kicked out :
# Pdav Psp62 Psp79 -> not available (yet?)
# Pant -> not published


# ## help : print this help
# rule help :
# 	shell :
# 		"sed -n 's/^##//p' Snakefile"

## all : all download all the genomes
rule all :
	input :
		expand("data/{sp}/genome.fa.gz", sp=species.split(' '))

#rule download_reads
#	input :
#		expand("data/{sp}/genome.fa.gz", sp=species.split(' '))

rule genome_download :
	output :
		"data/{sp}/genome.fa.gz"
	shell :
		"scripts/download.sh {output}"

# rule reads_downlaod :
# 	output :
# 		"data/{sp}/{SRA_id}_1.fastq.gz", "data/{sp}/{SRA_id}_2.fastq.gz"
# 	shell :
# 		"fastq-dump --accession {SRA_id} --outdir data/{sp} --split-files --gzip"
