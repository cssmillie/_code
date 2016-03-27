library(Seurat)
library(methods)

seur = readRDS('human.ming500.minc25.G125.MinPts25.seur.rds')

i = commandArgs(trailingOnly=T)[[1]]

print(paste(i, 'pca'))
cells.use = names(seur@ident[seur@ident == i])
seur = subsetData(seur, cells.use=cells.use)
seur = mean.var.plot(seur, y.cutoff=.1, y.high.cutoff=40, x.low.cutoff=.1, fxn.x=expMean, fxn.y=logVarDivMean, do.contour=FALSE)
seur = pca(seur, pcs.store=50)
print(paste(i, 'project'))
seur = jackStraw(seur, num.pc=50)
npcs = jackStraw.permutation.test(seur)$r
if(is.na(npcs) | npcs < 5){npcs = 5}
print(paste(i, 'npcs 1', npcs))
seur = project.pca(seur, pcs.store=50)
genes.use = pca.sig.genes(seur, (1:npcs), pval.cut=1e-3, max.per.pc=250)
print(paste(i, 'pca 2'))
seur = pca(seur, pc.genes=genes.use, pcs.store=50)
seur = jackStraw(seur, num.pc=50)
npcs = jackStraw.permutation.test(seur)$r
if(is.na(npcs) | npcs < 5){npcs = 5}
print(paste(i, 'npcs 1', npcs))
print(paste(i, 'tsne'))
seur = run_tsne(seur, dims.use=(1:npcs), do.fast=T, pt.size=.5)
saveRDS(seur, file=paste('human.ming500.minc25.G125.MinPts25.', i, '.seur.rds', sep=''))
