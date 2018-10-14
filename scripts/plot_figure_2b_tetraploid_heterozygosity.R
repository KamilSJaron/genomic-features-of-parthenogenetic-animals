# plotting barplots of heterozygosity of tetraploid speices
#
# Meloidogyne_javanica
# Meloidogyne_arenaria
# Adineta_vaga
# Adineta_ricciae
# Rotaria_macrura
# Rotaria_magnacalcarata

heter_tab <- read.table('tables/tetraploid_heterozygosity.tsv', header = T, row.names = 1)

heter_tab$heter <- (heter_tab$AABB + heter_tab$A)
heter_tab['Rmag1','heter'] <- heter_tab['Rmag1','AABB'] # all measured heterozygosity is between ohnologs
heter_tab['Rmag1','A'] <- 0 # truth is NA, close to 0, but don't know

ylim <- c(0, max(heter_tab$heter, na.rm = T))

heter_tab['Avag1','heter'] <- ceiling(max(heter_tab$heter, na.rm = T))

pal <- c("#D95F02", "#7570B3")

pdf('figures/fig2b_heterozygosity_of_tetraploids.pdf')

bar_pos <- barplot(heter_tab$heter, col = c('grey', rep(pal[2], 5)),
                   ylab = 'Heterozygosity [%]')
barplot(heter_tab$A, col = pal[1], add = T)

text(bar_pos, (heter_tab$heter - 1), c('high', heter_tab$AABB[-1]))
text(bar_pos, 0.25, c(heter_tab$A[1:3], '~0', heter_tab$A[5:6]))

sp_labels <- c('A. vaga', 'A. ricciae',
               'R. macrura', 'R. magnacalcarata',
               'M. javanica', 'M. arenaria')

sp_labels <- lapply(paste0(sp_labels, " "), function(x){bquote(italic(.(x)))})
text(bar_pos, -1, do.call(expression, sp_labels), las = 1, srt = 20, xpd = TRUE)

# legend('topright', bty = 'n', c('AABB', 'AA\' or AAAB'), pch = 20, col = rev(pal))

dev.off()