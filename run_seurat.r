library(optparse)

# Get input arguments
options = list(
	make_option('--dge', help='dge file', type='character'),
	make_option('--ming', help='number of genes per cell', type='integer'),
	make_option('--minc', help='number of cells per gene', type='integer'),
        make_option('--out', help='output prefix', type='character')
)
parser = OptionParser(option_list=options)
args = parse_args(parser);

library(methods)
library(Seurat)

# Make filenames
log_fn = paste(args$out, 'log.rds', sep='.')
seur_fn = paste(args$out, 'seur.rds', sep='.')

# Get log
log = list()
if(file.exists(log_fn)){log = readRDS(log_fn)}

# Check log
check_log = function(query){
    if(is.null(log[[query]])){
	return(FALSE)
    } else {
        return(TRUE)
    }
}

# Print status
print_status = function(line){
    print(paste(args$out, ": ", line, sep=''))
}

if(!check_log('seurat')){

# Load DGE (rows = genes, cols = barcodes)
print('Loading DGE')
data = read.table(args$dge, sep='\t', header=T, row.names=1)
print('Transforming data')
data = apply(data, 2, function(a){10000*a/sum(a)})
data = data.frame(log(data + 1))
print('Seurat object')
seur = new('seurat', raw.data = data)
seur = setup(seur, project=args$out, min.cells=args$minc, min.genes=args$ming, calc.noise=F, is.expr=0, names.delim='_', names.field=1)
print(seur)
print('Variable genes')
pdf(paste(args$out, 'mean_var.pdf', sep='.'))
seur = mean.var.plot(seur, y.cutoff=1, y.high.cutoff=40, x.low.cutoff=.1, fxn.x=expMean, fxn.y=logVarDivMean, do.contour=FALSE)
dev.off()
print('Save Seurat object')
saveRDS(seur, file=seur_fn)
log$seurat = 1
} else {
seur = readRDS(seur_fn)
}

# Run PCA + Jackstraw
if(!check_log('pca')){
print('PCA')
seur = pca(seur, pcs.store=50)
pdf(paste(args$out, 'pca.pdf', sep='.'))
pca.plot(seur, 1, 2, pt.size=1)
dev.off()
print('Jackstraw')
seur = jackStraw(seur, num.pc=50)
pc_sig = jackStraw.permutation.test(seur)
log$pc_sig = min(pc_sig, 50)
print('Save Seurat object')
saveRDS(seur, file=seur_fn)
log$pca = 1
saveRDS(log, file=log_fn)
}

# Calculate TSNE
if(!check_log('tsne')){
print('TSNE')
seur = run_tsne(seur, dims.use = 1:log$pc_sig$r, do.fast=T, pt.size=1)
pdf(paste(args$out, 'tsne.pdf', sep='.'), width=14, height=7)
tsne.plot(seur, pt.size=1, label.cex.text=.25)
dev.off()
print('Save Seurat object')
saveRDS(seur, file=seur_fn)
log$tsne = 1
saveRDS(log, file=log_fn)
}

# Projected PCA + Jackstraw
if(!check_log('pca2')){
print('Projected PCA')
seur = project.pca(seur, do.print=T, pcs.store=50)
genes.use = pca.sig.genes(seur, 1:log$pc_sig$r, pval.cut=1e-3, max.per.pc=250)
seur = pca(seur, pc.genes=genes.use, do.print=T, pcs.store=50)
pdf(paste(args$out, 'pca2.pdf', sep='.'))
pca.plot(seur, 1, 2, pt.size=1)
dev.off()
print('Projected Jackstraw')
seur = jackStraw(seur, num.pc=50)
pc_sig = jackStraw.permutation.test(seur)
log$pc_sig = min(pc_sig, 50)
print('Save Seurat object')
saveRDS(seur, file=seur_fn)
log$pca2 = 1
saveRDS(log, file=log_fn)
}

# Calculate TSNE
if(!check_log('tsne2')){
print('TSNE')
seur = run_tsne(seur, dims.use = 1:log$pc_sig$r, do.fast=T, pt.size=1)
pdf(paste(args$out, 'tsne2.pdf', sep='.'), width=14, height=7)
tsne.plot(seur, pt.size=1, label.cex.text=.25)
dev.off()
print('Save Seurat object')
saveRDS(seur, seur_fn)
log$tsne2 = 1
saveRDS(log, file=log_fn)
}

saveRDS(log, file=log_fn)
