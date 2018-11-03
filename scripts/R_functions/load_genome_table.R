load_genome_table <- function(columns = NA){

    genome_tab <- read.table('tables/genome_table.tsv',
                         header = T, stringsAsFactors = F, skip = 1, check.names = F)
    rownames(genome_tab) <- genome_tab$code
    # take only relevant columns and reverse the order of rows
    if ( ! any(is.na(columns)) ){
        genome_tab <- genome_tab[,columns]
    }

    genome_tab[is.na(genome_tab$hybrid_origin),'hybrid_origin'] <- "unknown"
    genome_tab[is.na(genome_tab$reproduction_mode),'reproduction_mode'] <- "unknown"

    hyb_origins <- c("no", "unknown", "yes")
    repr_modes  <- c("gamete_duplication", "terminal_fusion", "central_fusion", "unknown_automixis",
                     "unknown", "functional_apomixis")

    # ordering
    genome_tab$hybrid_origin <- ordered(genome_tab$hybrid_origin, levels=hyb_origins )
    genome_tab$reproduction_mode <- ordered(genome_tab$reproduction_mode, levels=repr_modes )
    genome_tab
}