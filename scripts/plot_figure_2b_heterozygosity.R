library('gsheet')
library('RColorBrewer')

############
# Get data #
############

literature_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1T4BHQxzMGMlWiNJ7G9OLPXzdNBCqMw0mSO7C_ggGEbc/edit?usp=sharing", format='csv'),
                            stringsAsFactors = F, skip = 1, header = T, check.names = F)
reported_heterozygosity <- literature_data[,c("code","haplotype_divergence[%]")]

genome_tab <- read.table('tables/genome_table.tsv',
                         header = T, stringsAsFactors = F, skip = 1, check.names = F)
rownames(genome_tab) <- genome_tab$code

genome_tab['Rmac1','heterozygosity'] <- 0.123 # these estimates should be added to the big genome table, not here
genome_tab['Rmag1','heterozygosity'] <- 0.492
genome_tab['Dcor1','reproduction_mode'] <- "unknown_automixis"

genome_tab$reproduction_mode[is.na(genome_tab$reproduction_mode)] <- "unknown"
genome_tab$hybrid_origin[is.na(genome_tab$hybrid_origin)] <- "unknown"
genome_tab$hybrid_origin[genome_tab$hybrid_origin == "possibly"] <- "suggested"
genome_tab$ploidy[is.na(genome_tab$ploidy)] <- "unknown"

hyb_origins <- c("no", "unknown", "suggested", "yes")
repr_modes  <- c("gamete_duplication", "terminal_fusion", "central_fusion", "unknown_automixis",
                 "apomixis", "unknown")

# ordering
genome_tab$hybrid_origin <- ordered(genome_tab$hybrid_origin, levels=hyb_origins )
genome_tab$reproduction_mode <- ordered(genome_tab$reproduction_mode, levels=repr_modes )

########
# plot #
########

pal <- c(rgb(0.9411, 0.8941, 0.2588), # yellow - gamete duplication
         rgb(0.9019, 0.6235, 0     ), # orange - terminal fusion
         rgb(0.8352, 0.3686, 0     ), # vermillion - central fusion
         rgb(0.8,    0.4745, 0.6549), # pink - automixis unknown
         rgb(0,      0.4470, 0.6980), # blue - apomixis
         rgb(0,      3/5,    2/5     )) # black - unknown

plot_ploits <- function(hyb_origin = "no", repr_mode = "gamete_duplication"){
         which_rm <- which(repr_mode == repr_modes)
         at <- which(hyb_origin == hyb_origins) +  which_rm / 7 - 0.5
         subset <- genome_tab[genome_tab$hybrid_origin == hyb_origin & genome_tab$reproduction_mode == repr_mode, c("heterozygosity", "ploidy")]
         if( nrow(subset) > 0 ){
                  points(rep(at, nrow(subset)), subset$heterozygosity, col = pal[which_rm], pch = ifelse(subset$ploidy == 2, 19, 15), cex = 1.3)
                  points(rep(at, nrow(subset)), subset$heterozygosity, pch = ifelse(subset$ploidy == 2, 21, 22), cex = 1.3)
         }
}

pdf('figures/fig2_heterozygosity.pdf', width=8, height=8)
# png('figures/fig2_heterozygosity.png')

# 'c(bottom, left, top, right)'
par(mar = c(4.5, 4.5, 2, 1) + 0.1)

with(genome_tab,
         boxplot(heterozygosity ~ hybrid_origin, bty = 'n', axes=FALSE,
                  xlab = "Hybrid origin", ylab = "Heterozygosity [%]",
                  pch = NA, col = rgb(0.878,0.878,0.878), cex.lab = 1.4
         )
)
axis(2, las=1)
axis(1, labels = hyb_origins, at = 1:4, tick = F, line = F, cex.axis = 1.4)

for ( ho in hyb_origins ){
         for ( rm in repr_modes ){
                  plot_ploits(ho, rm)
         }
}

legend('topleft', bty = 'n', c(repr_modes,'diplod','polyploid'), col = c(pal,1,1), pch = c(rep(20,6), 15), cex = 1.2)


dev.off()