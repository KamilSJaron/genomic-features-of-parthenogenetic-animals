# the naive expectation is that
# AB - AB * BC + BC = est_heterozygosity (which is the measured heterozygosty)
# B is the reference, A and C are alt genome copies and AB * BC is A * C distance, the random overlap
# AB (1 - BC) + BC = est_heterozygosity
# AB = (est_heterozygosity - BC) / (1 - BC)
## AB = BC (the most extreme case - equdistant)
# AB (2 - AB) = est_heterozygosity

get_exp_triallelic <- function(est_heterozygosity){
    # 2*x - x^2 = est_heterozygosity
    # x^2 - 2*x - est_heterozygosity = 0
    100 * ((1 - sqrt(1 - (est_heterozygosity / 100)))^2)
}

triploids <- read.table('tables/triploid_heterozygosity.tsv', header = T, stringsAsFactors=F)
observed_het <- (triploids$AAB + triploids$ABC)
naively_expected_trial <- sapply(observed_het, get_exp_triallelic)

# plot(triploids$ABC ~ naively_expected_ABC, xlim = c(0,1), ylim = c(0,1), pch = 20)
# text(naively_expected_ABC, triploids$ABC, triploids$code, pos = 1)
# lines(c(0,1),c(0,1))

total_heterozygosities <- seq(0, 8, len = 100)
triallelic_exp <- sapply(total_heterozygosities, get_exp_triallelic)

tiff("figures/Supp_fig2_rep_TE_patterns.tiff")

plot(triallelic_exp ~ total_heterozygosities,
     xlab = 'Heterozygosity [%]',
     ylab = 'Naively expected triallelic loci [%]', ylim = c(0, max(triploids$ABC)),
     type = 'l')
# points(observed_het, naively_expected_trial, pch = 20)
points(observed_het, triploids$ABC, pch = 20)
text(observed_het, triploids$ABC, triploids$code, pos = 1)

legend('topleft', bty = 'n', pch = c(NA, 20), lty = c(1, NA), c('expected', 'observed'))

dev.off()

######
# SIMULATION VERIFICATION
######

# genome_size = 10000
#
# generate_derived_haplotype <- function(ancestor, genome_size, mutaion_p){
#     muts <- rpois(1, round(mutaion_p * genome_size))
#     mut_locations <- round(runif(muts, 1, genome_size))
#     derived_seq <- ancestor
#     for (pos in mut_locations){
#         derived_seq[pos] <- sample(c('A', 'C', 'T', 'G'), 1)
#     }
#     derived_seq
# }
#
# calculate_divergence <- function(mutaion_p){
#     ancestor <- sample(c('A', 'C', 'T', 'G'), size = genome_size, replace = T)
#
#     h1 <- generate_derived_haplotype(ancestor, genome_size, mutaion_p)
#     h2 <- generate_derived_haplotype(ancestor, genome_size, mutaion_p)
#     h3 <- generate_derived_haplotype(ancestor, genome_size, mutaion_p)
#
#     AB_dis = sum(h1 != h2)
#     BC_dis = sum(h2 != h3)
#     AC_dis = sum(h1 != h3)
#
#     triallelic = sum(h1 != h2 & h2 != h3 & h1 != h3)
#     het = sum(((h1 != h2) + (h2 != h3) + (h1 != h3)) > 1)
#
#     c(mutaion_p, c(het, triallelic, AB_dis, BC_dis, AC_dis) / genome_size)
# }
#
# make_replicates <- function(mutaion_p, replicates){
#     results_mat <- sapply(rep(mutaion_p, replicates), calculate_divergence)
#     apply(results_mat, 1, median)
# }
#
# muts_to_simulate <- c(seq(0.001,0.01, by = 0.002),
#                       seq(0.01,0.04, by = 0.005),
#                       seq(0.05,0.3, by = 0.05))
# summary_per_mutation_rate <- as.data.frame(t(sapply(, make_replicates, 200)))
# colnames(summary_per_mutation_rate) <- c('mutation_p', 'het', 'triallelic', 'AB_dis', 'BC_dis', 'AC_dis')
#
# more_data <- as.data.frame(t(sapply(seq(0.01,0.04, by = 0.01), make_replicates, 200)))
# colnames(more_data) <- c('mutation_p', 'het', 'triallelic', 'AB_dis', 'BC_dis', 'AC_dis')
#
# summary_per_mutation_rate <- rbind(more_data, summary_per_mutation_rate)
#
# more_data <- as.data.frame(t(sapply(, make_replicates, 200)))
# colnames(more_data) <- c('mutation_p', 'het', 'triallelic', 'AB_dis', 'BC_dis', 'AC_dis')
#
#
# plot(c(sapply(heter, get_half)^2) ~ heter, type = 'l', col = 'green', xlab = 'heterozygosity [%]', ylab = 'naively expected overlap [%]', xlim = c(0, 0.1))
# points(summary_per_mutation_rate$het, summary_per_mutation_rate$triallelic, pch = 20)
