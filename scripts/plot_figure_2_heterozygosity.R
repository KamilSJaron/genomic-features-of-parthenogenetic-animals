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

# fix wrongly reported ploidy using using kmer spectra analysis
genome_tab[c('Mflo1', 'Mflo2'),'ploidy'] <- 3
# add unoknown ploidy using kmer spectra analysis
genome_tab[c('Ment1'),'ploidy'] <- 3

genome_tab$reproduction_mode[is.na(genome_tab$reproduction_mode)] <- "unknown"
genome_tab$hybrid_origin[is.na(genome_tab$hybrid_origin)] <- "unknown"
genome_tab$hybrid_origin[genome_tab$hybrid_origin == "possibly"] <- "suggested"
genome_tab$ploidy[is.na(genome_tab$ploidy)] <- "unknown"

hyb_origins <- c("no", "unknown", "suggested", "yes")
repr_modes  <- c("gamete_duplication", "automixis", "apomixis", "unknown")

# ordering
genome_tab$hybrid_origin <- ordered(genome_tab$hybrid_origin, levels=hyb_origins )
genome_tab$reproduction_mode <- ordered(genome_tab$reproduction_mode, levels=repr_modes )

########
# plot #
########

pal <- c(rgb(1,      1,      0.4274), # yellow - gamete duplication
         rgb(0,      0.5725, 0.5725), # drak green - automixis
         rgb(0,      0.4274, 0.8588), # blue - apomixis
         rgb(0.5725, 0,      0     )) # red - unknown

plot_ploits <- function(hyb_origin = "no", repr_mode = "gamete_duplication"){
         which_rm <- which(repr_mode == repr_modes)
         at <- which(hyb_origin == hyb_origins) +  which_rm / 5 - 0.5
         het <- genome_tab[genome_tab$hybrid_origin == hyb_origin & genome_tab$reproduction_mode == repr_mode, "heterozygosity"]
         if( length(het) > 0 ){
                  points(rep(at, length(het)), het, col = pal[which_rm], pch = 20, cex = 1.3)
                  points(rep(at, length(het)), het, pch = 1, cex = 1.3)
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

legend('topright', bty = 'n', repr_modes, col = pal, pch = 20, cex = 1.2)


dev.off()