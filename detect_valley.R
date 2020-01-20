#!/usr/bin/env Rscript
## Rscript to run DipSeek on the roadmap data
## To get help :
# Rscript --vanilla detect_valley.R --help

##########
# Parser #
##########

library("optparse")
library("data.table")
library("zoo")
library("tidyverse")

option_list = list(
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="valleys.out", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-l", "--read_length"), type="integer", 
              help="Average read extended length used for coverage computing", metavar="integer"),
  make_option(c("-n", "--testrun"), type="logical", default="FALSE", 
              help="If set as TRUE, will sample 10 peaks for test run [default= %default]", metavar="logical"),
  make_option(c("-c", "--chromosome"), type="character", default="*", 
              help="Regex of the chromosome you want to keep for the analysis", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}


#############
# Functions #
#############

# Make a function to split bed region larger than 30kb into smaller 10kb slidding windows
split_peaks <- function(df_in){
  large_peaks <- df_in %>% group_by(peak) %>% summarise(start = min(start), end = max(end)) %>% filter(end - start >= 30000) %>% select(peak)
  if (length(large_peaks$peak) == 0){return(df_in)
  } else {
    large_df <- c()
    for (p in seq(1,length(large_peaks$peak))){
      df <- df_in %>% filter(peak == large_peaks$peak[p])
      df$split <- (df$end - min(df$start)) %/% 10000
      df <- df %>% mutate(peak = paste0(peak,'.',split))
      df$split <- NULL
      large_df <- rbind(df, large_df)
    }
    small_df <- df_in %>% filter(!peak %in% large_peaks$peak)
    return(rbind(small_df, large_df))
    }
}

# Make a function to expand compress bedgraph format to single base pair format
expand_bedgraph <- function(df){
  df_out <- c()
  for (i in 1:nrow(df)){
    df_out <- rbind(df_out, data.frame('chromosome' = df[i,'chromosome'], 'bp' = seq(df[i,'start'],df[i,'end']-1),'cov' = df[i,'cov'],'peak'=df[i, 'peak']))
  }
  return(df_out)
}

# Make a functin to assess best smooting span value for a given peak
span_finder <- function(df_in){
  spanvalues = seq(0.05,1,0.05)
  mymetric <- matrix(ncol = 3, nrow = length(spanvalues))
  index = 1
  # Compute metrics for each span value
  for (i in spanvalues){
    theloess = loess(cov ~ bp, data = df_in, span = i)
    thecurve = predict(theloess)
    mymetric[index,] = c(i,theloess$s, sd(diff(thecurve)))
    index = index + 1
  }
  # Find the best span based on metrics values
  dfm <- data.frame("span" = mymetric[,1], "standard_error" = scale(mymetric[,2]), "roughness" = scale(mymetric[,3]))
  dist2d<-function(a,b,c){v1<- b - c; v2<- a - b; m<-rbind(v1,v2); d<-det(m)/sqrt(sum(v1*v1))}
  a2<-c(dfm$standard_error[dfm$standard_error==min(dfm$standard_error)][1],dfm$roughness[dfm$standard_error==min(dfm$standard_error)][1])
  b2<-c(dfm$standard_error[dfm$roughness==min(dfm$roughness)][1],dfm$roughness[dfm$roughness==min(dfm$roughness)][1])
  dfm$dist<-apply(dfm,1,function(x) dist2d(c(x[2],x[3]),a2,b2))
  best_span = dfm$span[dfm$dist==max(dfm$dist)]
  return(dfm[dfm$span == best_span,])
}


peakvalley_finder_updated <- function(df_in, best_span, length_read){
  # Run loess and find inflexion points and coverage sum 50bp around
  out = predict(loess(cov ~ bp, data = df_in, span = best_span),  seq(min(df_in$bp), max(df_in$bp), 1), se = TRUE)
  infl <- c(0, diff(diff(out$fit) > 0) , 0)
  landscape = data.frame(index=seq(min(df_in$bp), max(df_in$bp), 1), TF=infl, cov = df_in$cov, covaround = rollsum(df_in$cov, k = 100, fill = NA, align = 'center')/length_read) %>% filter(TF != 0)
  # Check if we have a valley-peak-valley pattern before we look at lag and lead
  if(all(landscape$TF == rep(landscape$TF[c(1,2)],length.out = length(landscape$TF)))){
    landscape <- landscape %>% mutate(lead_covaround = lead(covaround),
                                      lag_covaround = lag(covaround),
                                      lead_cov = lead(cov),
                                      lag_cov = lag(cov),
                                      local_max_dist = lead(index) - lag(index)) %>% filter(TF == 1 & !is.na(covaround))
    return(landscape)} else {warning(paste0("Peak ",df_in[1,'peak'],"does not have the expected P-V-P pattern"))}
}


#######
# Run #
#######

# Load data coverage and peak list
df_init <- fread(opt$file, header=FALSE, col.names = c('chromosome','start','end','cov','score','peak'), colClasses = c('character','integer','integer','integer','numeric','character'))
df_init <- split_peaks(df_init)
peaks_init <- unique(as.data.frame(df_init[,c('chromosome','peak')])) %>% filter(grepl(opt$chromosome, chromosome))
if (opt$testrun){peaks_init <- sample_n(peaks_init, 10)}

# Run bestspan and valley finder on each peak per chromosome
for (j in unique(peaks_init$chromosome)){ 
  print(paste0("processing ",j))
  peaks <- peaks_init$peak[peaks_init$chromosome == j]
  final_df <- c()
  final_peak <- c()
  for (i in seq(1,length(peaks))){
    best_span <- data.frame("span" = NA, "standard_error" = NA, "roughness" = NA)
    mypeak = peaks[i]
    mydf <- as.data.frame(df_init[which(df_init$peak == mypeak),])
    df <- expand_bedgraph(mydf)
    if (length(df$cov) > 500){
      suppressWarnings(best_span <- span_finder(df))
      if (best_span$span %in% seq(0.05,1,0.05)){
        df_out <- peakvalley_finder_updated(df, best_span$span, opt$read_length)
        if (nrow(df_out) != 0) {
          df_out <- df_out %>% mutate(peak = mypeak)
          final_df <- rbind(final_df, df_out)
          final_peak <- rbind(final_peak, data.frame(peak = mypeak, chromosome = mydf[1,'chromosome'], span = best_span$span, stderr = best_span$standard_error, roughness = best_span$roughness))
        }
      } else {print(paste("peak",mypeak,"does not output a proper span value"))}
    }
  }
  write.table(final_df, file = paste0(opt$out,'.',j,'.valleys.out'), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  write.table(final_peak, file = paste0(opt$out,'.',j,'.peaks.out'), col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
}
