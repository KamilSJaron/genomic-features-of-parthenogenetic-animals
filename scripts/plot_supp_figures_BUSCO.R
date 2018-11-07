############
# Get data #
############

source('scripts/R_functions/load_genome_table.R')

genome_tab <- load_genome_table(c(1:13))
genome_tab <- genome_tab[nrow(genome_tab):1,]

# kick out those that have no TEs estimated
genome_tab <- genome_tab[!is.na(genome_tab$complete), ]

####
# Adding genome size esitmates form data that are available

fillGenomeTab <- function(genome_tab, sp){
    cols <- genome_tab$species == sp
    het <- genome_tab[cols,'haploid_length[M]']
    genome_tab[which(cols)[is.na(het)],'haploid_length[M]'] <- mean(het, na.rm = T)
    genome_tab
}

genome_tab <- fillGenomeTab(genome_tab, 'Meloidogyne_arenaria')
genome_tab <- fillGenomeTab(genome_tab, 'Meloidogyne_javanica')
genome_tab <- fillGenomeTab(genome_tab, 'Meloidogyne_incognita')
genome_tab <- genome_tab[genome_tab$code != 'Dpul1',]

########
# duplications vs assembly
########

ploidies <- genome_tab[,"ploidy"]
haplotypes_in_asm <- genome_tab[,'assembly_size[M]'] / genome_tab[,'haploid_length[M]']

# I want to see that none of the assemblies is more than ploidy * haploid genome size
# and indeed diploid assemblies are 80 - 200% of gnome size
# triploid asm are 100 - 270% of the haploid genome size
# and tetraploid are 170 - 400%
# boxplot(haplotypes_in_asm ~ ploidies, xlab = 'ploidy', ylab = 'asm size / haploid genome size')

# png('figures/Supp_figure_BUSCO_duplicates.png')
pdf('figures/Supp_figure_BUSCO_duplicates.pdf')
    plot(haplotypes_in_asm ~ genome_tab$duplicated, pch = genome_tab$ploidy + 13,
         xlab = '% of duplicated BUSCO genes', ylab = 'assembly length / haploid genome size')
    text(genome_tab$duplicated, haplotypes_in_asm, genome_tab$code, pos = 1)
    legend('topright', legend = 2:4, bty = 'n', pch = 15:17, title = 'ploidy')
dev.off()

# png('figures/Supp_figure_BUSCO_missing_vs_N50.png')
pdf('figures/Supp_figure_BUSCO_missing_vs_N50.pdf')
    plot(log10(genome_tab[,'N50[k]']) ~ genome_tab$missing, pch = genome_tab$ploidy + 13,
         xlab = '% of missing BUSCO genes', ylab = 'log10( N50 of genome assembly)')
    text(genome_tab$missing, log10(genome_tab[,'N50[k]']), genome_tab$code, pos = 1)
    legend('topright', legend = 2:4, bty = 'n', pch = 15:17, title = 'ploidy')
dev.off()

########
# plot #
########

spaces <- c(0.1, 0.1, 1, rep(0.1, 15),1 , rep(0.15, 5), 1, rep(0.15, 3), 1)
busco_pal <- c('#fed976', '#fd8d3c','grey')

png('figures/Supp_figure_BUSCO.png')
pdf('figures/Supp_figure_BUSCO.pdf', width = 10, height = 6)
# plot BUSCO
# 'c(bottom, left, top, right)'
par(mar = c(5, 4.5, 1, 2) + 0.1)

bar_pos <- barplot(rep(100, nrow(genome_tab)),
                   xlab = "BUSCO [ % ]", cex.lab = 1.3,
                   col = busco_pal[3], horiz = T, axes = F, names.arg = F, space = spaces)
axis(1)
axis(2, at = bar_pos, labels = genome_tab$code, las = 1)

barplot(genome_tab$complete + genome_tab$fragmented, space = spaces,
        col = busco_pal[1], horiz = T, add = T, axes = F, names.arg = F)

# barplot(genome_tab$complete + genome_tab$fragmented, space = spaces,
#         col = busco_pal[2], horiz = T, add = T, axes = F, names.arg = F)
# barplot(genome_tab$complete + genome_tab$fragmented - genome_tab$duplicated, space = spaces,
#         col = busco_pal[1], horiz = T, add = T, axes = F, names.arg = F)

lines(c(90,90), range(bar_pos) + c(-0.5,0.5), lty = 3, lwd = 2)

dev.off()