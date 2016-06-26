library(optparse)

# This function plots a TSNE map with '--field1', '--clust' and '--label' options
# It also plots clusters from a different Seurat object, requiring:
# '--seur2' (to get cell ids)
# '--clust' (to get clusters)
# '--find' and '--replace' (map names from seur2 -> seur1 with gsub)

option_list = list(make_option('--seur1', help='seurat file for plotting (rds)', default=''),
                   make_option('--field1', help='color by ident field', type='integer', default=0),
                   make_option('--clust1', help='color by cluster (rds)', default=''),
		   make_option('--label1', help='labels for clust1 (tab)', default=''),
                   make_option('--seur2', help='color by another seurat object (rds)', default=''),
                   make_option('--clust2', help='color by another seurat object (rds)', default=''),
                   make_option('--label', help='label clusters?', default=FALSE, action='store_true'),
                   make_option('--search', help='search regex (seur2 -> seur1)', default=''),
                   make_option('--replace', help='replace regex (seur2 -> seur1)', default=''),
                   make_option('--both', help='plot both tsne maps?', default=FALSE, action='store_true'),
                   make_option('--out', help='output prefix')
                   )

if(interactive()){
    args = list(); args$seur1 = 'mistrg_bm.1000.50.combat.seur.rds'; args$seur2 = '../../data/mistrg_bm.500.25.seur.rds'; args$out = 'qq';
    args$field1 = 0; args$clust1 = ''; args$clust2 = ''; args$label = FALSE; args$search = ''; args$replace = ''; args$both = FALSE;
} else{
    args = parse_args(OptionParser(option_list=option_list))
}

library(Seurat)
library(methods)

# read seurat object
seur1 = readRDS(args$seur1)

# set ident with field1
if(args$field1 > 0){
    ident = sapply(strsplit(colnames(seur1@data), '\\.'), '[', field1)
    seur1 = set.ident(seur1, ident.use=ident)
}

# set ident with clust1
if(args$clust1 != ''){
    
    ident = readRDS(args$clust1)$membership
    
    if(args$label1 != ''){
        label1 = read.table(args$label1, sep='\t', row.names=1)
	ident = label1[ident,1]
    }
    
    seur1 = set.ident(seur1, ident.use=ident)
}

# set ident with another seurat object
if(args$seur2 != ''){
    # read seurat object
    seur2 = readRDS(args$seur2)
    # fix cell names
    if(args$search != '' & args$replace != ''){colnames(seur2@ident) = gsub(args$search, args$replace, colnames(seur2@ident))}
    # intersect cell names
    names = intersect(colnames(seur1@data), colnames(seur2@data))
    if(length(names) == 0){stop('regex')}
    # get cell identities
    if(args$clust2 != ''){
        ident = readRDS(args$clust2)$membership
        seur2 = set.ident(seur2, ident.use=ident)
    }
    # map to seurat object
    new_ident = data.frame(row.names=colnames(seur1@data), stringsAsFactors=F)
    if(args$label == FALSE){new_ident$x = NA}
    if(args$label == TRUE){new_ident$x = 'na'}
    new_ident[names,'x'] = as.character(seur2@ident[names])
    new_ident = as.character(new_ident$x)
    # set identities
    seur1 = set.ident(seur1, ident.use=new_ident)
}

if(args$both == FALSE){
    # plot single tsne
    pdf(paste(args$out, '.tsne.pdf', sep=''), width=12, height=8)
    tsne.plot(seur1, pt.size=.5, do.label=args$label)
    dev.off()
} else {
    # plot both tsnes
    pdf(paste(args$out, '.both.tsne.pdf', sep=''), width=16, height=8)
    par(mfrow=c(1,2))
    tsne.plot(seur1, pt.size=.5, do.label=args$label)
    tsne.plot(seur2, pt.size=.5, do.label=args$label)
    dev.off()
}
