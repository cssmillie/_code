# This script merges multiple DGEs into a single DGE

library(plyr)
library(stringr)

# Get input arguments (1=filenames, 2=ids, 3=gzfile, 4=genes/cell cutoff, 5=pdf)
args = commandArgs(trailingOnly=T)
fns = readLines(args[[1]])
ids = readLines(args[[2]])
cutoff = as.integer(args[[3]])
pdf_fn = args[[4]]
out_fn = args[[5]]

# Get list of DGEs
dges = llply(fns, function(a){read.table(a, header=T, row.names=1, sep='\t')})

# Relabel colnames to prevent collisions
for(i in 1:length(dges)){colnames(dges[[i]]) = paste(ids[i], colnames(dges[[i]]), sep='_')}

merge_dges = function(dges){
    # merge a list of dges into a single dge
    rows = sort(unique(unlist(llply(dges, rownames)))) # get rownames across all dges
    cols = sort(unique(unlist(llply(dges, colnames)))) # get colnames across all dges
    dge = data.frame(matrix(0, nrow=length(rows), ncol=length(cols))) # initialize dge matrix
    rownames(dge) = rows
    colnames(dge) = cols
    for(i in 1:length(dges)){dge[rownames(dges[[i]]), colnames(dges[[i]])] = dges[[i]]} # add all dges
    return(dge)
}

dge = merge_dges(dges)

# Plot the number of genes per cell
pdf(pdf_fn)
par(mfrow=c(1,2))
gene_count = sort(colSums(dge > 0), decreasing=T)
plot(gene_count, 1:length(gene_count), type='l', xlim=c(0,1500), xlab='genes per cell (cutoff)', ylab='number of cells')
plot(gene_count, cumsum(gene_count), type='l', xlim=c(0,1500), xlab='genes per cell (cutoff)', ylab='total genes')
dev.off()

# Filter DGE by gene count
dge = dge[, colSums(dge > 0) >= cutoff]

# Write DGE to outfile
write.table(dge, out_fn, quote=F, sep='\t')
