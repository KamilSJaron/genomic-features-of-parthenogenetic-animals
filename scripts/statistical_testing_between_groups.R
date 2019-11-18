genome_tab <- read.table('tables/genome_table.tsv',
                         header = T, stringsAsFactors = F, skip = 1, check.names = F)
rownames(genome_tab) <- genome_tab$code

# filter unknow hybrid origin
# genome_tab <- genome_tab[!is.na(genome_tab$hybrid_origin),]

#####
# I personally think that statistics is not well suited to the data we have
# we have limited number of mesurements and we also use biological knowledge
# about the samples to describe pattern we see,
# which is hard to encapsulate in statistical framework.
# moreover our samples are species with various degrees of phylogenetic relativness
# they are far from indipendent samples
#
# that sayed, stats kind of support the claims we made
#  - 69% of veriability in heterozygosity is explained simply by hybrid origin
#  - reproduction modes are untestable (therefore insignificant in statistics vocabulary)
#  - mann-whiney test of heterozygosity vs hybrid origin is also significant (less assumptions, p-value = 0.000999)
#  - TEs or Repeats are not significant to hybrid origin (p-value for both = 0.3636)

### testing heterozygosity
het_tab <- genome_tab[!is.na(genome_tab$heterozygosity), c('hybrid_origin', 'callular_mechanism', 'ploidy', 'heterozygosity')]
het_tab <- het_tab[!is.na(het_tab$hybrid_origin),]

wilcox.test(heterozygosity ~ hybrid_origin, het_tab)

# 	Wilcoxon rank sum test
#
# data:  heterozygosity by hybrid_origin
# W = 0, p-value = 0.000999
# alternative hypothesis: true location shift is not equal to 0

### parametric

summary(lm(heterozygosity ~ callular_mechanism, data = het_tab))

# nothing will be significant if there will be as many categories as

summary(lm(heterozygosity ~ hybrid_origin, data = het_tab))

# Call:
# lm(formula = heterozygosity ~ hybrid_origin, data = het_tab)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -3.7856 -0.2457  0.0290  0.4818  2.9844
#
# Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)
# (Intercept)        0.2460     0.8038   0.306 0.764814
# hybrid_originyes   5.2696     1.0025   5.256 0.000202 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
# Residual standard error: 1.797 on 12 degrees of freedom
# Multiple R-squared:  0.6972,	Adjusted R-squared:  0.672
# F-statistic: 27.63 on 1 and 12 DF,  p-value: 0.0002023

het_tab$ploidy = het_tab$ploidy - 2
summary(lm(heterozygosity ~ hybrid_origin + ploidy, data = het_tab))

# both significant:
#                  Estimate Std. Error t value Pr(>|t|)
# (Intercept)        0.2460     0.5757   0.427  0.67739
# hybrid_originyes   3.4476     0.8850   3.895  0.00250 **
# ploidy             2.0497     0.5822   3.521  0.00479 **

# ploidy here is an effect of every extra haplotype above diploidy (tri -> + 2.05, tetra -> + 4.1)

### testing repeats
rep_tab <- genome_tab[!is.na(genome_tab$repeats), c('hybrid_origin', 'callular_mechanism', 'repeats', 'TEs')]

rep_tab <- rep_tab[!is.na(rep_tab$callular_mechanism) & !is.na(rep_tab$TEs),]

wilcox.test(TEs ~ hybrid_origin, rep_tab)

# 	Wilcoxon rank sum test
#
# data:  TEs by hybrid_origin
# W = 32, p-value = 1
# alternative hypothesis: true location shift is not equal to 0

rep_tab$callular_mechanism[rep_tab$callular_mechanism != "functional_apomixis"] <- "automixis"

# extracted directly fromt the output files to resolve the tie
rep_tab['Lcla1','TEs'] <- 12.538
rep_tab['Anan1','TEs'] <- 12.524
wilcox.test(TEs ~ callular_mechanism, rep_tab)
# 	Wilcoxon rank sum test
#
# data:  TEs by callular_mechanism
# W = 55, p-value = 0.8446
# alternative hypothesis: true location shift is not equal to 0
