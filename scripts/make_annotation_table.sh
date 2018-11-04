
# it's reeeeaaally naughty oneliner

for annotation in data/*/annotation.gff3.gz; do
   sp=$(echo $annotation | cut -f 2 -d /); \
   printf "$sp\t" >> tables/gene_annotations.tsv; \
   zcat $annotation | \
   awk 'BEGIN{genes = 0; transripts = 0; mRNA = 0} { if($3=="mRNA") mRNA++; if($3 == "gene") genes++; if($3 == "transcript") transripts++; } END{ print genes "\t" transripts "\t" mRNA }' >> tables/gene_annotations.tsv;
done
