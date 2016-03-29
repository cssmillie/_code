# This script merges multiple DGEs into a single DGE

library(optparse)
library(plyr)
library(stringr)

# Get input arguments
option_list = list(make_option('--fns', help='filenames'),
                   make_option('--ids', help='ids'),
                   make_option('--out', help='output prefix', default='test')
                   )
args = parse_args(OptionParser(option_list=option_list))

# Get list of DGEs
dges = llply(fns, function(a){read.table(a, header=T, row.names=1, sep='\t')})




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
