# plotting barplots of heterozygosity of tetraploid speices
#
# Meloidogyne_javanica
# Meloidogyne_arenaria
# Adineta_vaga
# Adineta_ricciae
# Rotaria_macrura
# Rotaria_magnacalcarata

heter_tab <- read.table('tables/tetraploid_heterozygosity.tsv', header = T, row.names = 1)

heter_tab$heter <- (heter_tab$AABB + heter_tab$AAAB + heter_tab$AB)
heter_tab['Rmag1','heter'] <- heter_tab['Rmag1','AABB'] # all measured heterozygosity is between ohnologs
# truth is NA, close to 0, but don't know
# I need to set up these values so the the plot cleates boxes of 0 size
heter_tab['Rmag1',c('AB', 'AAAB')] <- 0
# we also don'jave reliable estimate of AAAB and AABB for A. vaga.
heter_tab$single_A_het <- heter_tab$AAAB + heter_tab$AB
heter_tab$single_A_het[1] <- heter_tab$AB[1]
heter_tab$AAAB[1] <- 0

ylim <- c(0, max(heter_tab$heter, na.rm = T))
# also I know that A. vaga has undetected high heterozygosty.
# So I will round up the highest other species
# because we can assume it's more
heter_tab['Avag1','heter'] <- ceiling(max(heter_tab$heter, na.rm = T))

# ColorBrewer2 palette (dark2, 3 colours)
pal <- c("#D95F02", "#7570B3", "#1B9E77")

pdf('figures/fig2b_heterozygosity_of_tetraploids.pdf')

bar_pos <- barplot(heter_tab$heter, col = pal[2],
                   ylab = 'Heterozygosity [%]')
barplot(heter_tab$single_A_het, col = pal[1], add = T)
barplot(heter_tab$AAAB, col = pal[3], add = T)

text(bar_pos, (heter_tab$single_A_het + 0.25), c('high', heter_tab$AABB[-1]))
text(bar_pos[c(1, 6)], heter_tab$AAAB[c(1, 6)] + 0.25, c(heter_tab$AB[c(1, 6)]))
text(bar_pos[c(5, 6)], 0.25, heter_tab$AAAB[c(5, 6)])

sp_labels <- c('A. vaga', 'A. ricciae',
               'R. macrura', 'R. magnacalcarata',
               'M. javanica', 'M. arenaria')

sp_labels <- lapply(paste0(sp_labels, " "), function(x){bquote(italic(.(x)))})
text(bar_pos, -1, do.call(expression, sp_labels), las = 1, srt = 20, xpd = TRUE)

# legend('topright', bty = 'n', c('AABB', 'AA\' or AAAB'), pch = 20, col = rev(pal))

dev.off()