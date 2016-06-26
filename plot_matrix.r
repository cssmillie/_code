library(optparse)

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
    args$clust = ''
    args$out = 'test'
} else{
option_list = list(make_option('--seur', help='seurat file'),
                   make_option('--matrix', help='gene expression matrix', default=''),
		   make_option('--labels', help='calculate mean across labels', default=''),
                   make_option('--tsne', help='plot combined tsne', default=FALSE, action='store_true'),
                   make_option('--tsne_k', help='number of clusters to diplay in combined tsne', default=25, type='integer'),
                   make_option('--tsne_i', help='plot individual tsnes', default=FALSE, action='store_true'),
                   make_option('--hmap', help='hmap plot', default=FALSE, action='store_true'),
                   make_option('--clust', help='clusters file (rds)', default=''),
                   make_option('--out', help='output prefix', default='')
                   )
args = parse_args(OptionParser(option_list=option_list))
}

library(gplots)
library(Seurat)
library(methods)

# Load data
seur = readRDS(args$seur)
data = seur@scale.data

# Get reference genes
refs = read.table(args$matrix, sep='\t', header=T, row.names=1)

# Set cell identities
if(args$clust != ''){
    ident.use = readRDS(args$clust)$membership
    seur = set.ident(seur, ident.use=ident.use)
}

# Intersect gene names
genes.use = intersect(rownames(data), rownames(refs))
data = data[genes.use,]
refs = refs[genes.use,]

# Get scores (rows=single cells, cols=reference cells)
library(proxy)
corr = cor(data, refs, method='pearson')

# Plot individual results as TSNE (colors = module scores)
if(args$tsne_i){
    d = data.frame(x=seur@tsne.rot[,1], y=seur@tsne.rot[,2])
    for(cell in colnames(corr)){
        d$z = scale(as.numeric(corr[,cell]))
        d$z[d$z < 0] = 0
        ggplot(d, aes(x=x, y=y)) + geom_point(aes(colour=z)) + scale_colour_gradient(low='#eeeeee', high='red') + theme_minimal()
        ggsave(paste(args$out, cell, 'tsne.pdf', sep='.'), width=12, height=8)
    }
}

# Plot combined results as TSNE (colors = best modules)
if(args$tsne){
    
    # Initialize data
    d = data.frame(x=seur@tsne.rot[,1], y=seur@tsne.rot[,2])
    
    # Get best cluster for each cell
    if(args$labels == ''){
        k = min(args$tsne_k, nrow(corr))
        ref_clusters = kmeans(t(corr), k)$cluster
        ref_clusters_names = aggregate(names(ref_clusters), by=list(ref_clusters), sample, 1)[,-1]
        ref_clusters = as.character(sapply(ref_clusters, function(i){ref_clusters_names[i]}))
    } else {
        ref_clusters = read.table(args$labels, sep='\t', row.names=1, stringsAsFactors=F)
	ref_clusters = ref_clusters[colnames(refs),1]
    }
    ident = as.character(ref_clusters[apply(corr, 1, which.max)])
    #ident = as.character(apply(corr, 1, function(a){tapply(a, ref_clusters[}))
    #ident = apply(corr, 1, function(a){names(which.max(tapply(a, ref_clusters, mean)))})
    print(aggregate(apply(corr,1,max), list(ident), mean))
    # Plot results
    seur = set.ident(seur, ident.use=ident)
    tsne.plot(seur, pt.size=1)
    ggsave(paste(args$out, 'combine', 'tsne.pdf', sep='.'), width=12, height=8)
}

# Plot results as heatmap (colors = module scores)
if(args$hmap){
    corr = aggregate(corr, by=list(seur@ident), mean)
    rownames(corr) = as.character(corr[,1])
    corr = corr[,-1]
    pdf(paste(args$out, 'hmap.pdf', sep='.'), width=12, height=8)
    heatmap.2(as.matrix(corr), Rowv=TRUE, Colv=NULL, scale='row', trace='none', col=colorRampPalette(c('purple', 'black', 'gold'))(64))
    dev.off()
}
