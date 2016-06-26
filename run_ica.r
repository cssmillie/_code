library(fastICA)
library(optparse)
library(Seurat)
library(methods)

option_list = list(make_option('--seur', help='Seurat file'),
                   make_option('--out', help='Output prefix'))
args = parse_args(OptionParser(option_list=option_list))

seur = readRDS(args$seur)
seur = mean.var.plot(seur, y.cutoff=1, y.high.cutoff=40, x.low.cutoff=.1, fxn.x=expMean, fxn.y=logVarDivMean, do.contour=FALSE)
seur = ica(seur, ics.store=75)
saveRDS(seur, file=paste(args$out, '.seur.rds', sep=''))
