library(MCL)
library(Seurat)
library(methods)

# Read input arguments
args = commandArgs(trailingOnly=T)
seur = args[[1]]
do.scale = as.integer(args[[2]])
do.spearman = as.integer(args[[3]])
out = args[[4]]

# Get output prefix
out = paste(out, do.scale, do.spearman, sep='.')

# Read input data
seur = readRDS(seur)

# Get similarity matrix
if(do.scale == 0){data = seur@data[,colnames(seur@scale.data)]}
if(do.scale == 1){data = seur@scale.data}
if(do.scale == 2){data = t(seur@pca.rot)}
if(do.spearman == 0){sim = cor(data)}
if(do.spearman == 1){sim = cor(data, method='spearman')}

# Cluster data
q = mcl(sim, addLoops=T)

# Set identities
seur = set.ident(seur, ident.use=q$Cluster)

# Plot
pdf(paste(out, '.pdf', sep=''), width=14, height=7)
tsne.plot(seur)
dev.off()

# Write clusters
saveRDS(q, file=paste(out, '.mcl.rds', sep=''))

# Get marker genes
markers = find_all_markers(seur, test.use='roc')
write.table(markers, file=paste(out, '.markers.txt', sep=''), sep='\t', quote=F)
