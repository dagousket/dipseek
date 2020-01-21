#!/usr/bin/env Rscript
## Rscript to run statistical stat following DipSeek run on the roadmap data
## To get help :
# Rscript --vanilla call_significant.R --help

##########
# Parser #
##########

library("optparse")
library("data.table")
library("tidyverse")
library("hexbin")
library("gridExtra")

option_list = list(
  make_option(c("-v", "--valleys"), type="character", default=NULL, 
              help="input valley file", metavar="character"),
  make_option(c("-p", "--peaks"), type="character", default=NULL, 
              help="input peak file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="valleys_tested.out", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-d", "--diagnosis"), type="character", default="diagnosis.pdf", 
              help="output diagnosis plot name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$valleys)){
  print_help(opt_parser)
  stop("Input file not provided.", call.=FALSE)
}


#############
# Functions #
#############


#######
# Run #
#######

# Load valley and peak data
valley_init <- fread(opt$valleys, header=FALSE, col.names = c('index','TF','cov','covaround','lead_covaround','lag_covaround','lead_cov','lag_cov','local_max_dist','peak'))
peak_init <- fread(opt$peaks, header=FALSE, col.names = c('peak','chromosome','span','stderr','roughness'))
df_init <- valley_init %>% left_join(peak_init, by = 'peak') %>% filter(local_max_dist > 200)
pdf(file = opt$diagnosis, onefile = TRUE, width = 15 , height = 7)

# Run diagnosis for span fitting
gspan1 <- ggplot(peak_init, aes(x = roughness, y = stderr)) + geom_hex(binwidth = c(0.1, 0.1)) + theme_light() + coord_fixed(ratio = 1)
gspan2 <- ggplot(peak_init, aes(x = roughness, y = span)) + geom_hex(binwidth = c(0.1, 0.05)) + theme_light()
gspan3 <- ggplot(peak_init, aes(x = stderr, y = span)) + geom_hex(binwidth = c(0.1, 0.05)) + theme_light()
lay <- matrix(c(1,2,1,3), byrow = TRUE, nrow = 2)
grid.arrange(arrangeGrob(gspan1),arrangeGrob(gspan2), arrangeGrob(gspan3), layout_matrix = lay)

# Run poisson test and diagnosis
df_init <- df_init %>% drop_na() %>% mutate(lambda = (lead_covaround + lag_covaround)/2, x = round(covaround), log2fc = log2(lambda/x))
df_init <- df_init %>% mutate(pvalue = ppois(x, lambda = lambda, lower.tail = TRUE), evalue = (pvalue * dim(df_init)[1]))
gval1 <- ggplot(df_init, aes(x = log2fc, y = pvalue)) + geom_hex() + theme_light() + scale_y_reverse()
gval2 <- ggplot(df_init, aes(x = pvalue)) + geom_histogram(binwidth = 0.05, color = 'white') + theme_light() + coord_cartesian(xlim = c(0,1))
lay <- matrix(c(1,2), byrow = TRUE, nrow = 1)
grid.arrange(arrangeGrob(gval1),arrangeGrob(gval2), layout_matrix = lay)

# Output bed file of the retained valleys
df_out <- df_init %>% mutate(start = index - 100, end = index + 100, score = abs(1000 - pvalue * 1000)) %>% select(chromosome, start, end, peak, score)
write.table(df_out, file = opt$out, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
dev.off()