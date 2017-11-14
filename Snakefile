# from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
# NCBI = NCBIRemoteProvider(email="kamiljaron+ncbi@gmail.com") # email required by NCBI to prevent abuse

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

rule download :
	input :
		"tables/download_table.tsv"
	output :
		"data/{sp}/genome.fa.gz"
	shell :
		"scripts/download.sh {output}"
