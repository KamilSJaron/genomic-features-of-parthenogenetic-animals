# from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
# NCBI = NCBIRemoteProvider(email="kamiljaron+ncbi@gmail.com") # email required by NCBI to prevent abuse
# import numpy

species_with_genomes='Cbir Avag Fcan Lcla Dcor Dpac Minc1 Minc2 Mjav Mare Mflo Dpul'
species_with_genomes='Cbir Avag Fcan Lcla Dcor Minc3 Mjav2 Mare2 Mflo Mflo2 Ment'
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
		expand("data/{sp}/genome.fa.gz", sp=species_with_genomes.split(' '))
		expand("data/{sp}/reads_R1.fq.gz", sp=species_with_genomes.split(' '))
		expand("data/{sp}/reads_R2.fq.gz", sp=species_with_genomes.split(' '))

rule download_genome :
	output :
		"data/{sp}/genome.fa.gz"
	shell :
		"scripts/download_genome.sh {output}"

rule downlaod_reads :
 	output :
 		"data/{sp}/reads_R1.fq.gz", "data/{sp}/reads_R2.fq.gz"
 	shell :
		"scripts/download_reads.sh {sp}"
