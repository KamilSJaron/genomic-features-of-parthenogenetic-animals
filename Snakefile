# from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
# NCBI = NCBIRemoteProvider(email="kamiljaron+ncbi@gmail.com") # email required by NCBI to prevent abuse
# import numpy

# Get from the data table all samples with
# genomes to species_with_genomes
# reads   to species_with_reads
species_with_genomes = []
species_with_reads = []

with open('tables/download_table.tsv') as tab :
	tab.readline()
	for textline in tab :
		line = textline.split()
		# print(line[1])
		# line[2] is location of genome
		if line[2] != 'NA' :
			# print(line[2] + 'is not NA, adding :' + line[0])
			species_with_genomes.append(line[0])
		# line[3] is SRA accession number
		if line[3] != 'NA' :
			# print(line[3] + 'is not NA, adding :' + line[0])
			species_with_reads.append(line[0])

# all_samples are unique values in array of merged samples with reads and genome
all_samples = list(set(species_with_genomes + species_with_reads))
# all_species is just a vector of unique entries of 4 letter substrings
all_species = list(set(map(lambda x: x[0:4], all_samples)))

mapping_files = []
theta_files = []
genome_stat_files = expand("data/{sp}/genome.stats", sp=species_with_genomes)
busco_files = expand("data/{sp}/busco", sp=species_with_genomes)
MUMmer_aln_files = expand("data/{sp}/MUMmer", sp=species_with_genomes)
genomescope_files = expand("data/{sp}/genomescope", sp=species_with_reads)

wind_size = 1000000
# we need to find all combinations of sequencing reads and references, so
# we iterate though all species
for spec in all_species :
	# iterate though all references of that species
	for ref in [x for x in species_with_genomes if x.startswith(spec)]:
		# iterate though all samples with reads given species
		for samp in [x for x in species_with_reads if x.startswith(spec)]:
			# for heterozygosity estimates
			theta_file = 'data/' + samp + '/' + ref + '_w' + str(wind_size) + '_theta_estimates.txt'
			theta_files.append(theta_file)
			# for
			mapping_file = 'data/' + samp + '/map_to_' + ref + '.bam'
			mapping_files.append(mapping_file)

### if environmental variable USE_LOCAL contains anything, it will compute on /scratch/local
cluster_script = os.environ.get("USE_LOCAL")
if cluster_script == None :
	cluster_script = ""
else :
	cluster_script = "scripts/use_local.sh "

## all
rule all :
	input : theta_files, genome_stat_files, busco_files

## calculate_busco
rule calculate_busco :
	input : busco_files

## calculate_selfalignments
rule calculate_selfalignment :
	input : MUMmer_aln_files

## calculate_heterozygosity_using_kmers
rule calculate_heterozygosity_using_kmers :
	input : genomescope_files

## calculate_genome_stats : calculate genome length, N50 and number of contigs of all genomes
rule calculate_genome_stats :
	input : genome_stat_files

## calculate_thetas : calculate theta estimates
rule calculate_thetas :
	input : theta_files

## map_all : map all samples reads to all reference genomes of the same species
rule map_all :
	input : mapping_files

## download_all : download all genomes and reads usnig information form tables/download_table.tsv
rule download_all :
	input :
		expand("data/{sp}/genome.fa.gz", sp=species_with_genomes),
		expand("data/{sp}/reads_R1.fq.gz", sp=species_with_reads)

## trimm_all : trimm all reads
rule trimm_all :
	input :
		expand("data/{sp}/reads-trimmed-pair1.fastq.gz", sp=species_with_reads)

## annotate_all_repeats : annotate repreats using reads and assembly size as a proxy for genome size; needs to run on dee-serv04
rule annotate_all_repeats :
	input :
		expand("data/{sp}/dnaPipeTE", sp=species_with_reads)

##
## help : print this help
rule help :
	shell :
		"sed -n 's/^##//p' Snakefile"

rule download_genome :
	threads : 1
	resources : mem=2000000, tmp=3000
	output : "data/{sp}/genome.fa.gz"
	shell : cluster_script + "scripts/download_data.sh {wildcards.sp} genome tables/download_table.tsv {output}"

rule download_proteins :
	threads : 1
	resources : mem=2000000, tmp=3000
	output : "data/{sp}/proteins.fa.gz"
	shell : cluster_script + "scripts/download_data.sh {wildcards.sp} proteins tables/download_table.tsv {output}"

rule download_annotation :
	threads : 1
	resources : mem=2000000, tmp=3000
	output : "data/{sp}/annotation.gff3.gz"
	shell : cluster_script + "scripts/download_data.sh {wildcards.sp} annotation tables/download_table.tsv {output}"

rule download_reads :
	threads : 1
	resources : mem=2000000, tmp=30000
	output : "data/{sp}/reads_R1.fq.gz"
	shell : cluster_script + "scripts/download_reads.sh {wildcards.sp} tables/download_table.tsv data/{wildcards.sp}/reads_R"

rule trim_reads :
	threads : 8
	resources : mem=80000000, tmp=150000
	input : "data/{sp}/reads_R1.fq.gz"
	output : "data/{sp}/reads-trimmed-pair1.fastq.gz"
	shell : cluster_script + "scripts/trim_reads.sh data/{wildcards.sp}/reads_R[1,2].fq.gz data/{wildcards.sp}/reads-trimmed"

rule index_reference :
	threads : 1
	resources : mem=20000000, tmp=10000
	input : "data/{reference}/genome.fa.gz"
	output : "data/{reference}/genome.fa.gz.bwt"
	shell :
		cluster_script + "scripts/index_genome.sh {input} data/{wildcards.reference}/genome.fa.gz."

rule map_reads :
	threads : 16
	resources : mem=104857600, tmp=40000
	input : "data/{reference}/genome.fa.gz.bwt", "data/{sample}/reads-trimmed-pair1.fastq.gz"
	output : "data/{sample}/map_to_{reference}.bam"
	shell :
		cluster_script + "scripts/map_reads.sh {wildcards.sample} {wildcards.reference} data/{wildcards.sample}/reads-trimmed-pair[1,2].fastq.gz data/{wildcards.reference}/genome.fa.gz.* {output}"

rule index_bam :
	threads : 1
	resources : mem=20000000, tmp=10000
	input : "data/{sample}/map_to_{reference}.bam"
	output : "data/{sample}/map_to_{reference}.bam.bai"
	shell : cluster_script + "scripts/index_bam.sh {input} {output}"

rule annotate_repeats :
	threads : 12
	resources : mem=150000000, tmp=30000
	input : "data/{sample}/reads-trimmed-pair1.fastq.gz"
	output : "data/{sample}/dnaPipeTE"
	shell : "scripts/annotate_repeats.sh {input} {wildcards.sample} {output}"

rule estimate_theta :
	threads : 1
	resources : mem=50000000, tmp=50000
	input : "data/{sample}/map_to_{reference}.bam", "data/{sample}/map_to_{reference}.bam.bai"
	output : "data/{sample}/{reference}_w{window_size}_theta_estimates.txt"
	shell :
		cluster_script + "scripts/est_theta.sh {wildcards.sample} {wildcards.reference} {wildcards.window_size} {input} {output}"

rule plot_all :
	threads : 1
	resources : mem=500000, tmp=5000
	input : expand(theta_files)
	output : "figures/species_heterozygosity.png"
	shell : "Rscript scripts/parse_thetas.R"

rule run_busco :
	threads : 16
	resources : mem=32000000, tmp=50000
	input : "data/{sp}/genome.fa.gz", "data/busco_ref/metazoa_odb9"
	output : "data/{sp}/busco"
	shell : cluster_script + "scripts/busco.sh {input} {output}"

rule get_busco_reference :
	threads : 1
	resources : mem=1000000, tmp=10000
	output : "data/busco_ref/metazoa_odb9"
	shell : """
		mkdir -p data/busco_ref && cd data/busco_ref
		wget http://busco.ezlab.org/datasets/metazoa_odb9.tar.gz
		tar -zxf metazoa_odb9.tar.gz
		rm metazoa_odb9.tar.gz
	"""

rule genome_stats :
	threads : 1
	resources : mem=1000000, tmp=10000
	input : "data/{sp}/genome.fa.gz"
	output : "data/{sp}/genome.stats"
	shell : "python3 scripts/fasta2genomic_stats.py {input} 1> {output}"

rule align_genome_to_itself :
	threads : 16
	resources : mem=100000000, tmp = 40000
	input : "data/{sample}/genome.fa.gz"
	output : "data/{sample}/MUMmer"
	shell : cluster_script + "scripts/MUMmer_selfaln.sh {input} {output}"

rule genome_profiling :
	threads : 16
	resources : mem=64000000, tmp = 60000
	input : "data/{sample}/reads-trimmed-pair1.fastq.gz"
	output : "data/{sample}/genomescope"
	shell : cluster_script + "scripts/GenomeScope.sh {input} data/{wildcards.sample}/reads-trimmed-pair2.fastq.gz {output}"
