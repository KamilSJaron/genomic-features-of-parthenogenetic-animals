# from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
# NCBI = NCBIRemoteProvider(email="kamiljaron+ncbi@gmail.com") # email required by NCBI to prevent abuse
# import numpy

# Get from the data table all samples with
# genomes to species_with_genomes
# reads   to species_with_reads
species_with_genomes = []
species_with_reads = []
species_with_annotation = []

raw_read_files = []
sample_accesions = dict()
busco_refs = dict()

with open('tables/download_table.tsv') as tab :
	tab.readline()
	for textline in tab :
		line = textline.split()
		sp = line[0]
		# print(line[1])
		# line[2] is location of genome
		if line[2] != 'NA' :
			# print(line[2] + 'is not NA, adding :' + line[0])
			species_with_genomes.append(sp)
		# line[3] is SRA accession number
		if line[3] != 'NA' :
			# print(line[3] + 'is not NA, adding :' + line[0])
			species_with_reads.append(sp)
			for lib in line[3].split(',') :
				raw_lib_file = 'data/' + sp + '/raw_reads/' + lib + '_1.fastq.gz'
				trimmed_lib_file = 'data/' + sp + '/trimmed_reads/' + lib + '-trimmed-pair1.fastq.gz'
				raw_read_files.append(raw_lib_file)
				sample_accesions[sp] = sample_accesions.get(sp, []) + [trimmed_lib_file]

species_with_reads_and_genomes = [sp for sp in species_with_reads if sp in species_with_genomes]

for sp in ["Dcor1", "Dpac1", "Minc1", "Minc2", "Minc3", "Mjav1", "Mjav2", "Mare1", "Mare2", "Mflo1", "Mflo2", "Ment1", "Pdav1", "Ps591", "Psam1"]:
	busco_refs[sp] = 'data/busco_ref/nematoda_odb9'

# all_samples are unique values in array of merged samples with reads and genome
all_samples = list(set(species_with_genomes + species_with_reads))
# all_species is just a vector of unique entries of 4 letter substrings
all_species = list(set(map(lambda x: x[0:4], all_samples)))

genome_stat_files = expand("data/{sp}/genome.stats", sp=species_with_genomes)
busco_files = expand("data/{sp}/busco", sp=species_with_genomes)

mapping_files = []
theta_files = []
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

localrules : help, all, calculate_busco, calculate_selfalignment, calculate_heterozygosity_using_kmers, calculate_genome_kmer_content, calculate_genome_stats, calculate_thetas, map_all, download_all, trimm_all

## all
rule all :
	input : theta_files, genome_stat_files, busco_files

## calculate_busco
rule calculate_busco :
	input : busco_files

## calculate_selfalignments
rule calculate_selfalignment :
	input : expand("data/{sp}/MUMmer", sp=species_with_genomes)

## calculate_heterozygosity_using_kmers
rule calculate_heterozygosity_using_kmers :
	input : expand("data/{sp}/genomescope", sp=species_with_reads)

## calculate_kmer_profiles_in_genome
rule calculate_genome_kmer_content :
	input : expand("data/{sp}/KAT", sp=species_with_reads_and_genomes)

## calculate_colinearity
rule calculate_colinearity :
	input : expand("data/{sp}/MCScanX", sp=species_with_genomes)

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
		raw_read_files

## trimm_all : trimm all reads
rule trimm_all :
	input :
		sample_accesions.values()

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

rule download_annotation :
	threads : 1
	resources : mem=2000000, tmp=3000
	output : "data/{sp}/annotation.gff3.gz"
	shell : cluster_script + "scripts/download_data.sh {wildcards.sp} annotation tables/download_table.tsv {output}"

rule download_reads :
	threads : 1
	resources : mem=2000000, tmp=30000
	output : "data/{sp}/raw_reads/{accesion}_1.fastq.gz"
	shell : cluster_script + "scripts/download_reads.sh {wildcards.sp} {wildcards.accesion} data/{wildcards.sp}/raw_reads/{wildcards.accesion}"

rule trim_reads :
	threads : 8
	resources : mem=80000000, tmp=150000
	input : "data/{sp}/raw_reads/{accesion}_1.fastq.gz"
	output : "data/{sp}/trimmed_reads/{accesion}-trimmed-pair1.fastq.gz"
	shell : cluster_script + "scripts/trim_reads.sh data/{wildcards.sp}/raw_reads/{wildcards.accesion}_[1,2].fastq.gz data/{wildcards.sp}/trimmed_reads/{wildcards.accesion}"

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
	input : "data/{reference}/genome.fa.gz.bwt", lambda wildcards: sample_accesions[wildcards.sample]
	output : "data/{sample}/map_to_{reference}.bam"
	shell :
		cluster_script + "scripts/map_reads.sh {wildcards.sample} {wildcards.reference} data/{wildcards.sample}/trimmed_reads/*.fastq.gz data/{wildcards.reference}/genome.fa.gz.* {output}"

rule index_bam :
	threads : 1
	resources : mem=20000000, tmp=10000
	input : "data/{sample}/map_to_{reference}.bam"
	output : "data/{sample}/map_to_{reference}.bam.bai"
	shell : cluster_script + "scripts/index_bam.sh {input} {output}"

rule annotate_repeats :
	threads : 12
	resources : mem=150000000, tmp=30000
	input : lambda wildcards: sample_accesions[wildcards.sample], "data/{sample}/genomescope"
	output : "data/{sample}/dnaPipeTE"
	shell : "scripts/annotate_repeats.sh data/{wildcards.sample}/trimmed_reads {wildcards.sample} {output}"

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
	shell : "script scripts/parse_thetas.R"

rule run_busco :
	threads : 16
	resources : mem=32000000, tmp=50000
	input : "data/{sp}/genome.fa.gz", lambda wildcards: busco_refs.get(wildcards.sp, "data/busco_ref/metazoa_odb9")
	output : "data/{sp}/busco"
	shell : cluster_script + "scripts/busco.sh {input} {output}"

rule get_busco_reference :
	threads : 1
	resources : mem=1000000, tmp=10000
	output : "data/busco_ref/{reference}_odb9"
	shell : """
		mkdir -p data/busco_ref && cd data/busco_ref
		wget http://busco.ezlab.org/datasets/{wildcards.reference}_odb9.tar.gz
		tar -zxf {wildcards.reference}_odb9.tar.gz
		rm {wildcards.reference}_odb9.tar.gz
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
	threads : 1
	resources : mem=4000000, tmp = 60000
	input : "data/{sample}/jellyfish"
	output : "data/{sample}/genomescope"
	shell : "scripts/genomescope.sh data/{wildcards.sample}/jellyfish {output}"

rule jellyfish :
	threads : 16
	resources : mem=64000000, tmp = 60000
	input : lambda wildcards: sample_accesions[wildcards.sample]
	output : "data/{sample}/jellyfish"
	shell : cluster_script + "scripts/jellyfish.sh data/{wildcards.sample}/trimmed_reads tables/download_table.tsv {output}"

rule kmer_genome_content :
	threads : 8
	resources : mem=64000000, tmp = 60000
	input : lambda wildcards: sample_accesions[wildcards.sample], "data/{sample}/genome.fa.gz"
	output : "data/{sample}/KAT"
	shell : cluster_script + "scripts/KAT.sh {wildcards.sample} data/{wildcards.sample}/trimmed_reads data/{wildcards.sample}/genome.fa.gz {output}"

rule colinearity_analysis :
	threads : 16
	resources : mem=64000000, tmp=60000
	input : "data/{sample}/genome.fa.gz", "data/{sample}/annotation.gff3.gz"
	output : "data/{sample}/MCScanX"
	shell : cluster_script + "scripts/MCScanX_palindromes.sh {wildcards.sample} {input} {output}"
