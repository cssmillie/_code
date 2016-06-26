source('~/aviv/code/var_genes.r')

args = commandArgs(trailingOnly=T)
counts = args[[1]]
prefix = args[[2]]
counts = read.table(counts, sep='\t', header=T, row.names=1)
meanCVfit(counts, prefix=prefix)
