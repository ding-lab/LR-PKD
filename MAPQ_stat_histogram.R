#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
print(args)

MAPQ_stat <- function(input_file) {
  
  MAPQ <- read.table(input_file, header = FALSE)
  MAPQ <- as.numeric(unlist(MAPQ))
  print("The variance of MAPQs is:")
  print(var(MAPQ))
  print("The sumary statistics of MAPQs is:")
  print(summary(MAPQ))

  pdf(file = paste(input_file,"histogram","pdf",sep = "."), width = 15,height = 12)
  # Histogram
  hist(MAPQ, freq = FALSE, col = "gray", xlab = "MAPQ", main = "The distribution of read MAPQs in PKD1")
  lines(density(MAPQ), col="blue", lwd=2) # add a density estimate with defaults
  dev.off()

}
  
MAPQ_stat(args[1])
