library(Seurat)
library(methods)
library(optparse)
library(gplots)

# Read input arguments
if(interactive()){
    args = list()
    args$seur = 'human.ming500.minc25.G125.MinPts25.seur.rds'
    args$matrix = '~/aviv/db/dmap/output/dmap_expression_mapped.tsv'
    args$type = 'sum'
    args$tsne = TRUE
    args$hmap = TRUE
    args$ident = ''
    args$out = 'test'
} else{
option_list = list(make_option('--seur', help='seurat file'),
                   make_option('--matrix', help='gene expression matrix', default=''),
                   make_option('--tsne', help='tsne plot', default=FALSE, action='store_true'),
                   make_option('--hmap', help='hmap plot', default=FALSE, action='store_true'),
                   make_option('--ident', help='cluster list (1=cell, 2=group)', default=''),
                   make_option('--out', help='output prefix', default='')
                   )
args = parse_args(OptionParser(option_list=option_list))
}

# Load data
seur = readRDS(args$seur)
data = seur@scale.data

# Get reference genes
refs = read.table(args$matrix, sep='\t', header=T, row.names=1)

# Set cell identities
if(args$ident != ''){
    ident.use = read.table(args$ident, sep='\t', row.names=1)
    ident.use = ident[colnames(data),1]
    seur = set.ident(seur, ident.use=ident.use)
}

# Intersect gene names
genes.use = intersect(rownames(data), rownames(refs))
data = data[genes.use,]
refs = refs[genes.use,]

# Get scores (rows=single cells, cols=reference cells)
corr = cor(data, refs)

# Plot results as TSNE (colors = module scores)
if(args$tsne){
    d = data.frame(x=seur@tsne.rot[,1], y=seur@tsne.rot[,2])
    for(cell in colnames(corr)){
        d$z = scale(as.numeric(corr[,cell]))
        d$z[d$z < 0] = 0
        ggplot(d, aes(x=x, y=y)) + geom_point(aes(colour=z)) + scale_colour_gradient(low='#eeeeee', high='red') + theme_minimal()
        ggsave(paste(args$out, cell, 'tsne.pdf', sep='.'), width=14, height=7)
    }
}

# Plot results as heatmap (colors = module scores)
if(args$hmap){
    corr = aggregate(t(corr), by=list(seur@ident), mean)
    rownames(corr) = as.character(corr[,1])
    corr = corr[,-1]
    pdf(paste(args$out, 'hmap.pdf', sep='.'), width=14, height=7)
    heatmap.2(as.matrix(corr), Rowv=TRUE, Colv=NULL, scale='row', trace='none', col=colorRampPalette(c('purple', 'black', 'gold'))(64))
    dev.off()
}
