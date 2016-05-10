library(optparse)

option_list = list(make_option('--seur', help='seurat file'),
                   make_option('--ident', help='cluster rds'),
                   make_option('--out', help='output prefix')
    )
args = parse_args(OptionParser(option_list=option_list))

if(file.exists(paste(args$out, '.seur.rds', sep=''))){stop('Seurat file exists')}

library(methods)
library(Seurat)

seur = readRDS(args$seur)
ident = sapply(strsplit(colnames(seur@data), '\\.'), '[', 2)
seur = set.ident(seur, ident.use=ident)
pdf(paste(args$out, '.tsne.data.pdf', sep=''), width=12, height=8)
tsne.plot(seur, pt.size=.75)
dev.off()
ident = readRDS(args$ident)$membership
seur = set.ident(seur, ident.use=ident)
pdf(paste(args$out, '.tsne.clust.pdf', sep=''), width=12, height=8)
tsne.plot(seur, pt.size=.75, do.label=T)
dev.off()
saveRDS(seur, paste(args$out, '.seur.rds', sep=''))
