library(optparse)

option_list = list(make_option('--seur', help='Seurat file'),
	           make_option('--field', help='Batch field', type='integer', default=1),
		   make_option('--types', help='Cell types', default=''),
		   make_option('--out', help='Output prefix'))

args = parse_args(OptionParser(option_list=option_list))

library(methods)
library(Seurat)
library(sva)

print('Reading data')
seur = readRDS(args$seur)

print('Batch')
batch=data.frame(batch=sapply(strsplit(colnames(seur@data), '\\.'), '[', args$field))

print('Cell types')
if(args$types != ''){
    types = read.table(args$types, sep='\t', header=T, row.names=1, stringsAsFactors=F)
    types = as.character(types[colnames(seur@data),1])
} else{
    types = c()
}

print('Design matrix')
design = data.frame(batch=batch, types=types)
model = model.matrix(~ as.factor(types), data=design)

print('Combating')
pdf(paste(args$out, '.combat.pdf', sep=''))
if(length(types) == 0){
    data = ComBat(dat=seur@data, batch$batch, par.prior=T, prior.plots=T)
} else {
    data = ComBat(dat=seur@data, batch$batch, par.prior=T, prior.plots=T, mod=model)
}
dev.off()

print('Write Seurat')
seur@data = data
seur@scale.data = scale(data)
saveRDS(seur, file=paste(out, '.combat.seur.rds', sep=''))
