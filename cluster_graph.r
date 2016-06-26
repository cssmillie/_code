library(Seurat)
library(methods)
library(optparse)
library(igraph)
library(cccd)

# Read input arguments
if(interactive()){
    args = list()
    args$seur = 'human.ming500.minc25.seur.rds'
    args$log = 'human.ming500.minc25.log.rds'
    args$pca = TRUE
    args$knn = 10
    args$type = 'louvain'
    args$out = 'test.out'
    args$test = TRUE
} else{
option_list = list(make_option('--seur', help='seurat file'),
                   make_option('--log', help='log file'),
                   make_option('--data', help='use data?', default=FALSE, action='store_true'),
                   make_option('--pca', help='use pca?', default=FALSE, action='store_true'),
                   make_option('--nnk', help='nearest neighbor k', type='integer'),
                   make_option('--type', help='louvain, infomap', default='louvain'),
                   make_option('--out', help='output prefix'),
                   make_option('--test', help='test', default=FALSE, action='store_true')
               )
args = parse_args(OptionParser(option_list=option_list))
}

# Read seurat object
seur = readRDS(args$seur)
log = readRDS(args$log)

# Get data
if(args$data == TRUE){
    d = t(seur@data)
} else if (args$pca == TRUE){
    num_pcs = ncol(seur@pca.rot)
    num_pcs = min(num_pcs, log$pc_sig)
    d = seur@pca.rot[,1:num_pcs]
} else {
    stop('must specify --data or --pca')
}

# Test dataset
if(args$test == TRUE){
    d = d[sample(1:nrow(d), 100, replace=F),]
}

# Build k-nearest-neighbor graph
g = nng(d, k=args$nnk)

# Calculate jaccard similarity
s = similarity(g, method='jaccard')

# Build weighted graph
g = graph.adjacency(s, mode='undirected', weighted=T)

# Cluster (Louvain method)
if(args$type == 'louvain'){
    g = cluster_louvain(g)
}
if(args$type == 'infomap'){
    g = cluster_infomap(g)
}

# Plot TSNE
seur = set.ident(seur, ident.use=g$membership)
pdf(paste(args$out, '.tsne.pdf', sep=''), width=12, height=8)
tsne.plot(seur, pt.size=.5, do.label=T)
dev.off()

# Write output
saveRDS(g, file=paste(args$out, '.clust.rds', sep=''))
