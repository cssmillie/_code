library(gplots)
library(methods)
library(optparse)
library(Seurat)

# Get input arguments
options = list(
    make_option('--seur', help='seurat file'),
    make_option('--clust', help='cluster rds (infomap or louvain)'),
    make_option('--plot', help='plot?', default=FALSE, action='store_true'),
    make_option('--out', help='output prefix')
)
parser = OptionParser(option_list=options)
args = parse_args(parser)

# Read data
seur = readRDS(args$seur)
clust = readRDS(args$clust)
seur = set.ident(seur, ident.use=clust$membership)
dmap = read.table('~/aviv/db/dmap/output/dmap_expression_mapped.tsv', sep='\t', header=T, row.names=1)

# Mean expression in each cluster
data = aggregate(t(seur@scale.data), list(seur@ident), mean)
data = t(data[,-1])

# Align data
genes.use = intersect(rownames(data), rownames(dmap))
data = data[genes.use,]
dmap = dmap[genes.use,]

# Correlation
corr = cor(data, dmap)

# Heatmap
if(args$plot == TRUE){
colsep = c(10, 14, 18, 27, 60, 72, 76, 89, 98, 103, 109, 119, 128, 148, 166, 211)
jpeg(paste(args$out, 'dmap.jpg', sep='.'), h=10, w=25)
heatmap.2(corr, Rowv=T, Colv=NULL, trace='none', colsep=colsep, col=colorRampPalette(c("purple","black","gold"))(64))
dev.off()
}

# Find marker genes - bimod
markers.bimod = find_all_markers(seur, test.use='bimod')
write.table(markers.bimod, file=paste(args$out, 'markers.bimod.txt', sep='.'), sep='\t', quote=F)
markers.use = subset(markers.bimod, avg_diff > 0)$gene
if(args$plot == TRUE){
jpeg(paste(args$out, 'bimod.jpg', sep='.'))
doHeatMap(seur, genes.use=markers.use, slim.col.label=T, remove.key=T, cexRow=.1)
dev.off()
}

# Find marker genes - roc
markers.roc = find_all_markers(seur, test.use='roc')
write.table(markers.roc, file=paste(args$out, 'markers.roc.txt', sep='.'), sep='\t', quote=F)
markers.use = subset(markers.roc, avg_diff > 0)$gene
if(args$plot == TRUE){
jpeg(paste(args$out, 'roc.jpg', sep='.'))
doHeatMap(seur, genes.use=markers.use, slim.col.label=T, remove.key=T, cexRow=.1)
dev.off()
}

# Find marker genes - tobit
markers.tobit = find_all_markers(seur, test.use='tobit')
write.table(markers.tobit, file=paste(args$out, 'markers.tobit.txt', sep='.'), sep='\t', quote=F)
markers.use = subset(markers.tobit, avg_diff > 0)$gene
if(args$plot == TRUE){
jpeg(paste(args$out, 'tobit.jpg', sep='.'))
doHeatMap(seur, genes.use=markers.use, slim.col.label=T, remove.key=T, cexRow=.1)
dev.off()
}
