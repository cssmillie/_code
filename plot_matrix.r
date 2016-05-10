library(Seurat)
library(methods)
library(optparse)
library(gplots)

# Read input arguments
if(interactive()){
    args = list()
    args$seur = '~/aviv/human/tsne/human.500.25.seur.rds'
    args$matrix = '~/aviv/db/dmap/output/dmap_expression_mapped.tsv'
    args$type = 'sum'
    args$tsne = TRUE
    args$tsne_k = 25
    args$tsne_i = FALSE
    args$hmap = FALSE
    args$ident = ''
    args$out = 'test'
} else{
option_list = list(make_option('--seur', help='seurat file'),
                   make_option('--matrix', help='gene expression matrix', default=''),
                   make_option('--tsne', help='plot combined tsne', default=FALSE, action='store_true'),
                   make_option('--tsne_k', help='number of clusters to diplay in combined tsne', default=25, type='integer'),
                   make_option('--tsne_i', help='plot individual tsnes', default=FALSE, action='store_true'),
                   make_option('--hmap', help='hmap plot', default=FALSE, action='store_true'),
                   make_option('--ident', help='clusters file (rds)', default=''),
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
    ident.use = readRDS(args$ident)$membership
    seur = set.ident(seur, ident.use=ident.use)
}

# Intersect gene names
genes.use = intersect(rownames(data), rownames(refs))
data = data[genes.use,]
refs = refs[genes.use,]

# Get scores (rows=single cells, cols=reference cells)
corr = cor(data, refs)

# Plot individual results as TSNE (colors = module scores)
if(args$tsne_i){
    d = data.frame(x=seur@tsne.rot[,1], y=seur@tsne.rot[,2])
    for(cell in colnames(corr)){
        d$z = scale(as.numeric(corr[,cell]))
        d$z[d$z < 0] = 0
        ggplot(d, aes(x=x, y=y)) + geom_point(aes(colour=z)) + scale_colour_gradient(low='#eeeeee', high='red') + theme_minimal()
        ggsave(paste(args$out, cell, 'tsne.pdf', sep='.'), width=14, height=7)
    }
}

# Plot combined results as TSNE (colors = best modules)
if(args$tsne){
    
    # Initialize data
    d = data.frame(x=seur@tsne.rot[,1], y=seur@tsne.rot[,2])
    k = min(args$tsne_k, ncol(corr))
    
    # Get best cluster for each cell
    ref_clusters = kmeans(t(corr), k)$cluster
    ref_clusters_names = aggregate(names(ref_clusters), by=list(ref_clusters), sample, 1)[,-1]
    ref_clusters = as.character(sapply(ref_clusters, function(i){ref_clusters_names[i]}))
    ident = as.character(ref_clusters[apply(corr, 1, which.max)])
    
    # Plot results
    seur = set.ident(seur, ident.use=ident)
    tsne.plot(seur, pt.size=1)
    ggsave(paste(args$out, 'combine', 'tsne.pdf', sep='.'), width=14, height=7)
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
