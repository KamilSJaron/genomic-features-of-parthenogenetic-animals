# # # # # # # #
# WILD SCRIPT #
# # # # # # # #

tripliods_gs <- read.table('tables/triploid_heterozygosity.tsv', header = T)

smudge_files <- paste0('data/', tripliods_gs$code, '/smudgeplot/', tripliods_gs$code, '_summary_table.tsv')

read.smudgesummary <- function(file){
    read.table(file, sep = '\t', skip = 1,
               col.names = c('smudge', 'kmer_pairs', 'proportion', 'center_sum', 'center_frac'))
}

triploid_smudges <- lapply(smudge_files, read.smudgesummary)

smudges_AAB <- sapply(triploid_smudges, function(x){ x[x$smudge == 'AAB', 'proportion'] })
smudges_AB <- sapply(triploid_smudges, function(x){ x[x$smudge == 'AB', 'proportion'] })
smudges_AB <- c(0,unlist(smudges_AB))

gs_div_score <- tripliods_gs$AB / (tripliods_gs$AB + tripliods_gs$AAB)
smudge_div_score <- smudges_AB / (smudges_AB + smudges_AAB)

png('figures/triploid_haplotype_structure.png')

par(mar = c(5,5,4,2))
plot(smudge_div_score ~ gs_div_score, pch = 20,
     xlab = 'GenomeScope \n ABC / (AAB + ABC)', ylab = 'smudgeplot \n AB / (AAB + ABC)')
text(gs_div_score, smudge_div_score, as.character(tripliods_gs$code), pos = 1)

title('Higher the value is, more equidistant the three haplotypes are')

dev.off()

# inter = 0.1
# intra = 0.005
#
# AABB = inter - (inter * intra * 2)
# AABC = (inter * intra * (1 - intra) * 2) + (1 - inter) * intra * intra
# AAAB = ((1 - inter) * intra * (1 - intra) * 2)
# ABCD = inter * intra * intra
#
# AABB + AABC + AAAB + ABCD
# AAAB
#
# inter = 0.06
# intra = 0.015
#
# AABB = inter - (inter * intra * 2)
# AABC = (inter * intra * (1 - intra) * 2) + (1 - inter) * intra * intra
# AAAB = ((1 - inter) * intra * (1 - intra) * 2)
# ABCD = inter * intra * intra
#
# AABB + AABC + AAAB + ABCD
# AAAB