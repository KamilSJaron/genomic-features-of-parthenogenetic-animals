# plotting barplots of heterozygosity of tetraploid speices
#
# Meloidogyne_javanica
# Meloidogyne_arenaria
# Adineta_vaga
# Adineta_ricciae
# Rotaria_macrura
# Rotaria_magnacalcarata

trip_tab <- read.table('tables/triploid_heterozygosity.tsv', header = T, row.names = 1)
tetra_tab <- read.table('tables/tetraploid_heterozygosity.tsv', header = T, row.names = 1)

trip_tab$heter <- (trip_tab$AAB + trip_tab$ABC)
tetra_tab$heter <- (tetra_tab$AABB + tetra_tab$AAAB + tetra_tab$AB)

# we also don'jave reliable estimate of AAAB and AABB for A. vaga.
tetra_tab$single_A_het <- tetra_tab$AAAB + tetra_tab$AB
tetra_tab$single_A_het[1:4] <- tetra_tab$AB[1:4]
tetra_tab$AAAB[1:4] <- 0

tetra_tab <- tetra_tab[c(5:6, 1:4),]

ylim <- c(0, max(tetra_tab$heter, na.rm = T))
# also I know that A. vaga has undetected high heterozygosty.
# So I will round up the highest other species
# because we can assume it's more
tetra_tab[3:6,'heter'] <- tetra_tab[3:6,'AABB']
# ceiling(max(tetra_tab$heter, na.rm = T)) + 1

# ColorBrewer2 palette (dark2, 3 colours)
# pal <- c("#D95F02", "#7570B3", "#1B9E77")
pal <- c("#D81B60CA", "#1E88E5", "#FFC107", "#004D40AA")

pdf('figures/fig3_heterozygosity_of_tetraploids.pdf', width = 8, height = 6)

spaces <- c(rep(0.4, 3), rep(0.2, 2), 0.6, 0.2, 0.4, rep(0.2, 3))
bar_pos <- barplot(c(rep(NA, 5), tetra_tab$heter), col = pal[1],
                   ylab = 'Decomposed heterozygosity [%]', space = spaces)
barplot(c(trip_tab$heter, tetra_tab$single_A_het), col = pal[3], add = T, space = spaces)
barplot(c(rep(NA, 5), tetra_tab$AB), col = pal[2], add = T, space = spaces)
barplot(c(trip_tab$ABC, rep(NA, 6)), col = pal[2], add = T, space = spaces)

shift <- 0.8
# 1 digit rounding
text(bar_pos, c(rep(NA, 5), tetra_tab$single_A_het + shift), c(rep(NA, 5), round(tetra_tab$AABB, 1)))
text(bar_pos[c(8, 9)], 0 + shift, round(c(tetra_tab$AB[c(3, 4)]), 1))
text(bar_pos[c(6, 7)], tetra_tab$AB[c(1, 2)] + shift, round(tetra_tab$AAAB[c(1, 2)], 1))
text(bar_pos[c(1:5)], trip_tab$ABC + shift, round(trip_tab$AAB, 1))
# text(bar_pos[c(2,3,5)], 0.25, round(trip_tab$ABC[c(2,3,5)], 1))


sp_labels <- c('P. virginalis', 'P. davidi',
               'M. incognita', 'M. floridensis', 'M. enterolobii',
               'M. javanica', 'M. arenaria',
               'A. vaga', 'A. ricciae',
               'R. macrura', 'R. magnacalcarata')

sp_labels <- lapply(paste0(sp_labels, " "), function(x){bquote(italic(.(x)))})
text(bar_pos, -2, do.call(expression, sp_labels), las = 1, srt = 20, xpd = TRUE)

# for (bar in c(1,2,5,7)) {
#     lty = ifelse(bar == 5, 1, 2)
#     lines(rep(mean(bar_pos[bar:(bar+1)]), 2), c(-5, 33), lty = lty)
# }


# legend('topright', bty = 'n', c('AABB', 'AA\' or AAAB'), pch = 20, col = rev(pal))

dev.off()