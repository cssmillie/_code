library(optparse)

# Read input arguments
option_list = list(make_option('--seur', help='seurat file'),
                   make_option('--log', help='log file'),
                   make_option('--out', help='output prefix', default='')
                   )
args = parse_args(OptionParser(option_list=option_list))

library(Seurat)
library(methods)
library(gplots)

# Load data
seur = readRDS(args$seur)

# Get PCs
pcs = seur@pca.rot

# Get number of genes
num_genes = colSums(seur@data > 0)

# Plot the number of genes
library(ggplot2)
d = data.frame(x=seur@tsne.rot[,1], y=seur@tsne.rot[,2], z=num_genes)
ggplot(d, aes(x=x,y=y)) + geom_point(aes(col=z)) + scale_colour_gradient(low='#eeeeee', high='red') + theme_minimal()
ggsave(paste(args$out, '.reads.tsne.pdf', sep=''), width=12, height=8)

# Find PCs correlated with the number of genes
for(j in 1:ncol(pcs)){
    print(cor(pcs[,j], num_genes))
}
