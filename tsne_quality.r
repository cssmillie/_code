library(optparse)

# Read input arguments
if(interactive()){
    args = list()
    args$seur = '/home/unix/csmillie/aviv/analysis/april2016_2/human/tsne/frozen.500.25.seur.rds'
    args$dge = '/home/unix/csmillie/aviv/analysis/april2016_2/human/dge/human.dge.txt.gz'
    args$out = 'test'
} else{
option_list = list(make_option('--seur', help='seurat file'),
                   make_option('--dge', help='dge file'),
                   make_option('--out', help='output prefix', default='')
                   )
args = parse_args(OptionParser(option_list=option_list))
}

library(Seurat)
library(methods)

# Load data
seur = readRDS(args$seur)
dge = read.table(args$dge, sep='\t', header=T, row.names=1)
dge = dge[,colnames(seur@data)]

# Set identity
ident = colnames(seur@data)
ident = gsub('\\.[^\\.]*$', '', ident)
seur = set.ident(seur, ident.use=ident)

# Quality
num_reads = colSums(dge)
num_genes = colSums(dge > 0)
quality = cbind(num_reads, num_genes)

# Plot TSNE
pdf(paste(args$out, '.collections.pdf', sep=''), width=12, height=8)
tsne.plot(seur, pt.size=1)
dev.off()
d = data.frame(x=seur@tsne.rot[,1], y=seur@tsne.rot[,2])
d$z = scale(num_reads)
ggplot(d, aes(x=x, y=y)) + geom_point(aes(colour=z)) + scale_colour_gradient(low='#eeeeee', high='red') + theme_minimal()
ggsave(paste(args$out, '.num_reads.pdf', sep=''), width=12, height=8)
d$z = scale(num_genes)
ggplot(d, aes(x=x, y=y)) + geom_point(aes(colour=z)) + scale_colour_gradient(low='#eeeeee', high='red') + theme_minimal()
ggsave(paste(args$out, '.num_genes.pdf', sep=''), width=12, height=8)
