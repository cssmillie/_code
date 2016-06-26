library(optparse)
option_list = list(make_option('--dge', help='dge file', default=''),
	           make_option('--seur', help='seurat file', default=''),
	      	   make_option('--matrix', help='expression matrix', default='/home/unix/csmillie/aviv/db/dmap/output/dmap_expression_mapped.tsv'),
		   make_option('--labels', help='expression labels', default='/home/unix/csmillie/aviv/code/cell_types/dmap.labels.txt'),
		   make_option('--mouse', help='mouse defaults', default=FALSE, action='store_true'),
		   make_option('--cutoff', help='correlation cutoff', type='numeric', default=0),
		   make_option('--out', help='output prefix')
		   )
args = parse_args(OptionParser(option_list=option_list))

library(Seurat)
library(methods)

# Load data
if(args$dge != ''){
    data = read.table(args$dge, sep='\t', header=T, row.names=1)
    data = apply(data, 2, function(a){10000*a/sum(a)})
    data = data.frame(log(data + 1))
    data = scale(data)
} else {
    seur = readRDS(args$seur)
    data = seur@scale.data
}

# Mouse data
if(args$mouse == TRUE){
    args$matrix = '/home/unix/csmillie/aviv/db/immgen/immgen.mm10.txt'
    args$labels = '/home/unix/csmillie/aviv/code/cell_types/immgen.labels.txt'
}

# Reference matrix
refs = read.table(args$matrix, sep='\t', header=T, row.names=1)

# Intersect gene names
genes.use = intersect(rownames(data), rownames(refs))
data = data[genes.use,]
refs = refs[genes.use,]

# Cluster scores
scores = cor(data, refs, method='pearson')

# Cluster labels
labels = read.table(args$labels, sep='\t', row.names=1, stringsAsFactors=F)
labels = labels[colnames(refs), 1]

# Map cells
ident = apply(scores, 1, function(a){
    q = tapply(a, labels, max)
    if(max(q) < args$cutoff){
        return(NA)
    } else{
        return(names(which.max(q)))
    }
})
ident = data.frame(ident=ident, row.names=colnames(data))

# Write results
write.table(ident, paste(args$out, '.cell_types.txt', sep=''), sep='\t', quote=F)
write.table(table(ident), paste(args$out, '.cell_types.summary.txt', sep=''), sep='\t', quote=F)

# Plot results
if(args$seur != ''){
    seur = set.ident(seur, ident.use=ident[colnames(seur@data),1])
    pdf(paste(args$out, '.cell_types.tsne.pdf', sep=''), width=12, height=8)
    tsne.plot(seur, pt.size=1)
    dev.off()
}
