#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


library(dplyr)
library(broom)

HybWW <- read.table(args[1], header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
HybWW2 <- HybWW %>% 
 group_by(NR) %>%
  do(tidy(binom.test(.$Counts_B73, .$Counts_B73+.$Counts_NAM, alternative = "two.sided")))
  
HybWW$p.value = HybWW2$p.value

HybWW$fdr =
      p.adjust(HybWW$p.value,
               method = "fdr")
               
write.table(HybWW, file = args[2], sep = "\t",
            row.names = TRUE, col.names = NA)

