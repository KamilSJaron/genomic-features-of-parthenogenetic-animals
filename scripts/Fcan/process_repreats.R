springtail_repeats <- read.table('data/Fcan/fcand_repeats.txt', header = T)

TEs_lower <- springtail_repeats$Repeat %in% c('SINE', 'LINE', 'LTR', 'DNA', 'Retroposon')
TEs_upper <- springtail_repeats$Repeat %in% c('SINE', 'LINE', 'LTR', 'DNA', 'Retroposon', 'Unknown')
overall <- 0.233

sum(springtail_repeats[TEs_lower,'Percentage']) * overall
sum(springtail_repeats[TEs_upper,'Percentage']) * overall

