library(optparse)

# Read input arguments
if(interactive()){
    args = list()
    args$seur = '~/aviv/human/tsne/human.500.25.seur.rds'
    args$pairs = 'pairs.in'
    args$out = 'test'
} else{
option_list = list(make_option('--seur', help='seurat file'),
                   make_option('--pairs', help='pairs list (1=gene1, 2=gene2)', default=''),
                   make_option('--out', help='output prefix', default='')
                   )
args = parse_args(OptionParser(option_list=option_list))
}

library(Seurat)
library(methods)
library(gplots)

# Load data
seur = readRDS(args$seur)
data = seur@data

# Get reference genes
pairs = read.table(args$pairs, stringsAsFactors=F)

# Intersect gene names
genes.use = intersect(rownames(data), unique(c(pairs[,1], pairs[,2])))
if(length(genes.use) == 0){stop('genes not found')}
data = data[genes.use,,drop=F]

# Get module scores (modules x cells)
d = data.frame(x=seur@tsne.rot[,1], y=seur@tsne.rot[,2])
for(i in 1:nrow(pairs)){
    gene1 = pairs[i,1]
    gene2 = pairs[i,2]
    if(gene1 %in% rownames(data) & gene2 %in% rownames(data)){
        d$z = as.numeric(data[gene1,]/(data[gene1,] + data[gene2,]))
        ggplot(d, aes(x=x, y=y)) + geom_point(aes(colour=z, alpha=.75)) + scale_colour_gradient2(low='blue', mid='purple', high='red', midpoint=.5) + theme_minimal()
        ggsave(paste(args$out, gene1, gene2, 'tsne.pdf', sep='.'), width=14, height=7)
    } else {print(paste('could not find genes', gene1, gene2))}
}