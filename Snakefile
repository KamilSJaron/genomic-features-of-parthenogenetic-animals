# from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
# NCBI = NCBIRemoteProvider(email="kamiljaron+ncbi@gmail.com") # email required by NCBI to prevent abuse
# import numpy

species_with_genomes='Cbir1 Avag1 Fcan1 Lcla1 Dcor1 Dpac1 Minc1 Minc2 Mjav1 Mare1 Mflo1 Dpul1 Hduj1'.split(' ')
species_with_reads='Cbir1 Avag1 Fcan1 Lcla1 Dcor1 Minc3 Mjav2 Mare2 Mflo1 Mflo2 Ment1 Hduj1'.split(' ')
# kicked out :
# Dpac -> not available reads in SRA
# Pdav Psp62 Psp79 -> not available (yet?)
# Pant -> not published
# to be added
# Dpul


all_samples = list(set(species_with_genomes + species_with_reads))
all_species = list(set(map(lambda x: x[0:4], all_samples)))

mapping_files = []
# iterate though all species
for species in all_species :
	# iterate though all references of that species
	for reference in list(filter(lambda x:species in x, species_with_genomes)):
		# iterate though all samples with reads given species
		for sample in list(filter(lambda x:species in x, species_with_reads)):
			mapping_file = 'data/' + sample + '/map_to_' + reference + '.bam'
			mapping_files.append(mapping_file)



# ## help : print this help
# rule help :
# 	shell :
# 		"sed -n 's/^##//p' Snakefile"

# rule map_all :
# 	input :


rule download_all :
	input :
		expand("data/{sp}/genome.fa.gz", sp=species_with_genomes),
		expand("data/{sp}/reads_R1.fq.gz", sp=species_with_reads)

rule download_genome :
	output : "data/{sp}/genome.fa.gz"
	shell : "scripts/download_genome.sh {output}"

rule downlaod_reads :
	output : "data/{sp}/reads_R1.fq.gz"
	shell : "scripts/download_reads.sh {output}"

rule index_reference :
	threads: 1
	resources: M=20000000
	input : "data/{reference}/genome.fa.gz"
	output : "data/{reference}/genome.fa.gz.bwt"
	shell :
		"scripts/cluster.sh scripts/index_genome.sh {input} data/{wildcards.reference}/genome.fa.gz."
#
# rule map_reads :
# 	input :  "data/{reference}/genome.fa.gz", "data/{reference}/genome.fa.bai", "data/{sample}/reads_R1.fq.gz"
# 	output : "data/{sample}/map_to_{reference}.bam"
