library(gplots)
library(methods)
library(optparse)
library(Seurat)

# Get input arguments
options = list(
	make_option('--out', help='output prefix', type='character'),
        make_option('--dbg', help='dbclust g parameter', type='numeric')
)
parser = OptionParser(option_list=options)
args = parse_args(parser)

# Read data
seur = readRDS(paste(args$out, '.seur.rds', sep=''))
log = readRDS(paste(args$out, '.log.rds', sep=''))
dmap = read.table('~/data/dmap/output/dmap_expression_mapped.tsv', sep='\t', header=T, row.names=1)

# Cluster
print('Cluster PCA')
seur = DBclust_dimension(seur, 1, 2, G.use=args$dbg, set.ident=T)

# Plot TSNE
pdf(paste(args$out, '.lab_tsne.pdf', sep=''))
tsne.plot(seur, pt.size=1, do.label=T)
dev.off()

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
colsep = c(10, 14, 18, 27, 60, 72, 76, 89, 98, 103, 109, 119, 128, 148, 166, 211)
pdf(paste(args$out, 'dmap.pdf', sep='.'), h=10, w=25)
heatmap.2(corr, Rowv=T, Colv=NULL, trace='none', colsep=colsep, col=colorRampPalette(c("purple","black","gold"))(64))
dev.off()

# Find marker genes - bimod
markers.bimod = find_all_markers(seur, test.use='bimod')
write.table(markers.bimod, file=paste(args$out, 'markers.bimod.txt', sep='.'), sep='\t', quote=F)
markers.use = subset(markers.bimod, avg_diff > 0)$gene
pdf(paste(args$out, 'bimod.pdf', sep='.'))
doHeatMap(seur, genes.use=markers.use, slim.col.label=T, remove.key=T, cexRow=.1)
dev.off()

# Find marker genes - roc
markers.roc = find_all_markers(seur, test.use='roc')
write.table(markers.roc, file=paste(args$out, 'markers.roc.txt', sep='.'), sep='\t', quote=F)
markers.use = subset(markers.roc, avg_diff > 0)$gene
pdf(paste(args$out, 'roc.pdf', sep='.'))
doHeatMap(seur, genes.use=markers.use, slim.col.label=T, remove.key=T, cexRow=.1)
dev.off()

# Find marker genes - tobit
markers.tobit = find_all_markers(seur, test.use='tobit')
write.table(markers.tobit, file=paste(args$out, 'markers.tobit.txt', sep='.'), sep='\t', quote=F)
markers.use = subset(markers.tobit, avg_diff > 0)$gene
pdf(paste(args$out, 'tobit.pdf', sep='.'))
doHeatMap(seur, genes.use=markers.use, slim.col.label=T, remove.key=T, cexRow=.1)
dev.off()
