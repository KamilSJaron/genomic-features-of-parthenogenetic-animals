
get_value <- function(line){
    as.numeric(strsplit(stat_lines[line], '\t')[[1]][[2]])
}

# function that unlists strsplit
ssplit <- function (s, split = "="){
    unlist(strsplit(s, split = split))
}

### BUSCO
read_busco <- function(busco_file){
    if( is.na(busco_file) ){
        return( rep(NA, 4) )
    }
    busco_file <- readLines(busco_file)
    total_genes <- as.numeric(ssplit(busco_file[15], '\t')[2])
    bscores <- c(complete = as.numeric(ssplit(busco_file[10], '\t')[2]),
                 fragmented = as.numeric(ssplit(busco_file[13], '\t')[2]),
                 duplicated = as.numeric(ssplit(busco_file[12], '\t')[2]),
                 missing = as.numeric(ssplit(busco_file[14], '\t')[2]))
    bscores <- round(100 * (bscores / total_genes), 2)
    return(bscores)
}

### GenomeScope
parse_genomescope_summary <- function(file){
    genoscope_file <- readLines(file)
    model_ploidy <- grepl("p = ", genoscope_file)
    est_ploidy <- as.numeric(ssplit(genoscope_file[model_ploidy], ' ')[3])

    homozygous_pattern <- paste(rep('A', est_ploidy), collapse='')
    line <- ssplit(genoscope_file[grepl(homozygous_pattern, genoscope_file)], ' ')
    if ( length(line) == 0 ){
        homozygous_pattern <- paste(rep('a', est_ploidy), collapse='')
        line <- ssplit(genoscope_file[grepl(homozygous_pattern, genoscope_file)], ' ')
    }
    line <- line[line != '']
    heterozygosity_position <- which(grepl(homozygous_pattern, line)) + 1
    heterozygosity_min <- round(100 - as.numeric(substr(line[heterozygosity_position], 0, nchar(line[heterozygosity_position]) - 1)), 2)
    heterozygosity_max <- round(100 - as.numeric(substr(line[heterozygosity_position + 1], 0, nchar(line[heterozygosity_position + 1]) - 1)), 2)
    heterozygosity <- mean(c(heterozygosity_min, heterozygosity_max))

    line <- ssplit(genoscope_file[grepl("Haploid", genoscope_file)], ' ')
    line <- gsub(",", "", line[line != ''])
    haploid_genome <- as.numeric(line[c(4, 6)])

    line <- ssplit(genoscope_file[grepl("Repeat", genoscope_file)], ' ')
    line <- gsub(",", "", line[line != ''])
    repeats <- round((mean(as.numeric(line[c(4, 6)])) * 100) / mean(haploid_genome), 2)

    return( c(round(mean(haploid_genome) / 1e6, 1),
              repeats,
              heterozygosity) )
}