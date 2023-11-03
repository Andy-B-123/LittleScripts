#!/apps/R/4.0.5/bin/Rscript

# Quick and dirty R script to generate coverage summaries from mosdepth thresholds output
# Default it's taking 1x coverage as the threshold
# First arguement is input filepath for mosdepth output file
# Second arguement is output filepath for processed data

library(data.table)
library(dplyr)
library(tidyverse)
library(plyr)
library(janitor)

args <- commandArgs(trailingOnly = TRUE)

#args_1 <- c("W:\\Harm_P450_paper_FINAL\\temp\\mosdepth\\29.Bedops.ByPangenomeCDS.thresholds.bed.gz")
#args_2 <- c("W:\\Harm_P450_paper_FINAL\\temp\\mosdepth\\29.Bedops.ByPangenomeCDS.regions.bed.gz")
#args_3 <- c("W:\\Harm_P450_paper_FINAL\\temp\\mosdepth\\29.outputBase")

sample_thresholds <- fread(args[1])

sample_thresholds
sample_thresholds_clean <- clean_names(sample_thresholds)
length(unique(sample_thresholds_clean$region))

sample_thresholds_clean$length <- sample_thresholds_clean$end - sample_thresholds_clean$start
sample_thresholds_clean_slim <- sample_thresholds_clean %>%
  select(c('number_chrom','start','end','region','length', everything())) %>%
  select(-c(start,end))

sample_thresholds_clean_grouped <- sample_thresholds_clean_slim %>%
  select(-c("number_chrom")) %>%
  group_by(region)

sample_thresholds_clean_grouped_summary <- sample_thresholds_clean_grouped %>%
  select(-region) %>%
  summarise_all(funs(sum))

sample_thresholds_clean_grouped_summary[, 3:ncol(sample_thresholds_clean_grouped_summary)] <- sapply(sample_thresholds_clean_grouped_summary[, 3:ncol(sample_thresholds_clean_grouped_summary)], function(x) x / sample_thresholds_clean_grouped_summary[, 2])

simple_summary <- sample_thresholds_clean_grouped_summary %>%
  select(region,length,x1x) %>%
  mutate(rounded_1x = round_any(x1x,0.1,f = round))

output_summary <- simple_summary %>%
  separate(region, into = c('ID',"parent"),sep = ';') %>%
  mutate(ID = gsub('ID=','',ID)) %>%
  mutate(parent = gsub("Parent=",'',parent)) %>%
  select(ID, parent, rounded_1x)

fwrite(output_summary,paste(args[3],".transCovExtent.tsv",sep = ""),na = 0,row.names = F, quote=FALSE )

### Below process the mosdepth output to get per transcript average coverage

sample_regions <- fread(args[2])

sample_regions
sample_regions_clean <- clean_names(sample_regions)
length(unique(sample_regions_clean$region))

sample_regions_clean$length <- sample_regions_clean$v3 - sample_regions_clean$v2


sample_thresholds_clean_grouped <- sample_regions_clean %>%
  select(-c(v1,v2,v3)) %>%
  select(v4,length,v5) %>%
  mutate(v5 = length * v5) %>%
  group_by(v4)

sample_thresholds_clean_grouped_summary <- sample_thresholds_clean_grouped %>%
  select(-v4) %>%
  summarise_all(funs(sum)) %>%
  mutate(trans_cov = v5 / length) %>%
  select(v4,trans_cov)

output_summary_cov <- sample_thresholds_clean_grouped_summary %>%
  separate(v4, into = c('ID',"parent"),sep = ';') %>%
  mutate(ID = gsub('ID=','',ID)) %>%
  mutate(parent = gsub("Parent=",'',parent)) %>%
  mutate(trans_cov = round(trans_cov, digits = 3)) %>%
  select(ID, parent, trans_cov)

fwrite(output_summary_cov,paste(args[3],".transCovAvg.tsv",sep = ""),na = 0,row.names = F, quote=FALSE )
