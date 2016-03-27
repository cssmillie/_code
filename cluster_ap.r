library(MCL)
library(Seurat)
library(methods)

# Read input arguments
args = commandArgs(trailingOnly=T)
seur = args[[1]]
sim = args[[2]]
out = args[[3]]

# Get output prefix
out = paste(out, do.scale, do.spearman, sep='.')

# Read input data
seur = readRDS(seur)

# Cluster data
if(sim == 0){s = negDistMat(r=1)}
if(sim == 1){s = corSimMat(r=1)}
if(sim == 2){s = linKernel(normalize=T)}
q = apcluster(s, t(seur@scale.data), 

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
