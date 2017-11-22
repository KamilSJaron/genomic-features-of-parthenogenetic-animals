# from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
# NCBI = NCBIRemoteProvider(email="kamiljaron+ncbi@gmail.com") # email required by NCBI to prevent abuse
# import numpy

species_with_genomes='Cbir1 Avag1 Fcan1 Lcla1 Dcor1 Dpac1 Minc1 Minc2 Mjav1 Mare1 Mflo1 Dpul1 Hduj1'
species_with_reads='Cbir1 Avag1 Fcan1 Lcla1 Dcor1 Minc3 Mjav2 Mare2 Mflo1 Mflo2 Ment1 Hduj1'
# kicked out :
# Dpac -> not available reads in SRA
# Pdav Psp62 Psp79 -> not available (yet?)
# Pant -> not published
# to be added
# Dpul

# ## help : print this help
# rule help :
# 	shell :
# 		"sed -n 's/^##//p' Snakefile"

## all : all download all the genomes
rule all :
	input :
		expand("data/{sp}/genome.fa.gz", sp=species_with_genomes.split(' ')), expand("data/{sp}/reads_R1.fq.gz", sp=species_with_reads.split(' '))

rule download_genome :
	output : "data/{sp}/genome.fa.gz"
	shell : "scripts/download_genome.sh {output}"

rule downlaod_reads :
 	output : "data/{sp}/reads_R1.fq.gz"
 	shell : "scripts/download_reads.sh {output}"
