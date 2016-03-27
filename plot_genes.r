library(optparse)

# Read input arguments
if(interactive()){
    args = list()
    args$seur = '~/aviv/bm/tsne/human.ming500.minc25.seur.rds'
    args$gene = 'SELP'
    args$modules = ''
    #args$modules = '~/aviv/db/markers/markers.test'
    args$type = 'sum'
    args$tsne = TRUE
    args$hmap = TRUE
    args$ident = ''
    args$out = 'test'
} else{
option_list = list(make_option('--seur', help='seurat file'),
                   make_option('--gene', help='gene name', default=''),
                   make_option('--modules', help='genes list (1=module, 2=gene, 3=direction)', default=''),
                   make_option('--type', help='sum, mean', default='sum'),
                   make_option('--tsne', help='tsne plot', default=FALSE, action='store_true'),
                   make_option('--hmap', help='hmap plot', default=FALSE, action='store_true'),
                   make_option('--ident', help='cluster list (1=cell, 2=group)', default=''),
                   make_option('--out', help='output prefix', default='')
                   )
args = parse_args(OptionParser(option_list=option_list))
}



library(Seurat)
library(methods)
library(gplots)

# Load data
seur = readRDS(args$seur)
data = seur@scale.data

# Get reference genes
if(args$gene != ''){
    refs = data.frame(cbind(args$gene, args$gene, 1), stringsAsFactors=F)
    refs[,3] = as.numeric(refs[,3])
}
if(args$modules != ''){
    refs = read.table(args$modules, sep='\t', stringsAsFactors=F)
}
if(ncol(refs) == 1){
    refs = as.data.frame(cbind(refs[,1], refs[,1], 1))
    refs[,3] = as.numeric(refs[,3])
}
if(ncol(refs) == 2){
    refs = as.data.frame(cbind(refs, 1))
}

# Set cell identities
if(args$ident != ''){
    ident.use = read.table(args$ident, sep='\t', row.names=1)
    ident.use = ident[colnames(data),1]
    seur = set.ident(seur, ident.use=ident.use)
}

# Intersect gene names
genes.use = intersect(rownames(data), unique(refs[,2]))
data = data[genes.use,,drop=F]
refs = refs[refs[,2] %in% genes.use,,drop=F]

# Get module scores (modules x cells)
scores = aggregate(1:nrow(refs), by=list(refs[,1]), function(i){
    q = data[refs[i,2],,drop=F] * refs[i,3]
    if(args$type == 'sum'){
        return(colSums(q))
    }
    if(args$type == 'mean'){
        return(colMeans(q))
    }
})
scores = as.data.frame(scores[[2]], row.names=scores[[1]])

# Plot results as TSNE (colors = module scores)
if(args$tsne){
    d = data.frame(x=seur@tsne.rot[,1], y=seur@tsne.rot[,2])
    for(module in rownames(scores)){
        d$z = scale(as.numeric(scores[module,]))
        d$z[d$z < 0] = 0
        ggplot(d, aes(x=x, y=y)) + geom_point(aes(colour=z)) + scale_colour_gradient(low='#eeeeee', high='red') + theme_minimal()
        ggsave(paste(args$out, module, 'tsne.pdf', sep='.'), width=14, height=7)
    }
}

# Plot results as heatmap (colors = module scores)
if(args$hmap){
    scores = aggregate(t(scores), by=list(seur@ident), mean)
    rownames(scores) = as.character(scores[,1])
    scores = scores[,-1]
    pdf(paste(args$out, 'hmap.pdf', sep='.'), width=14, height=7)
    heatmap.2(as.matrix(scores), Rowv=TRUE, Colv=NULL, trace='none', col=colorRampPalette(c('purple', 'black', 'gold'))(64))
    dev.off()
}
