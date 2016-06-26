library(optparse)

if(interactive()){
    args = list()
    args$seur1 = '~/aviv/mouse/tsne/mouse.1000.10.seur.rds'
    args$seur2 = '~/aviv/data/mouse1/cd45/CD45neg.1000.10.seur.rds'
    args$clust = '~/aviv/data/mouse1/cd45/CD45neg.1000.10.infomap.50.clust.rds'
    args$label = TRUE
    args$find = '^'
    args$replace = 'mouse1_cdn.'
    
} else{
option_list = list(make_option('--seur1', help='seurat file for plotting (rds)', default=''),
                   make_option('--seur2', help='seurat file for labeling (rds)', default=''),
                   make_option('--clust', help='cluster file (rds)', default=''),
                   make_option('--label', help='label clusters?', default=FALSE, action='store_true'),
                   make_option('--find', help='find pattern (for cluster names)', default=''),
                   make_option('--replace', help='replace pattern (for cluster names)', default='')
                   )
args = parse_args(OptionParser(option_list=option_list))
}

library(Seurat)
library(methods)

seur1 = readRDS(args$seur1)

if(args$seur2 == ''){
    if(args$clust != ''){
        ident = readRDS(args$clust)$membership
        seur1 = set.ident(seur1, ident.use = ident)
    }
}

if(args$seur2 != ''){
    if(args$clust == ''){stop('--seur2 needs --clust')}
    seur2 = readRDS(args$seur2)
    ident = readRDS(args$clust)$membership
    names = colnames(seur2@data)
    if(args$find != '' & args$replace != ''){
        names = gsub(args$find, args$replace, names)
    }
    if(sum(!names %in% colnames(seur1@data)) > 0){stop('regex')}
    new_ident = data.frame(row.names=colnames(seur1@data), stringsAsFactors=F)
    if(args$label == FALSE){new_ident$x = NA}
    if(args$label == TRUE){new_ident$x = 'na'}
    new_ident[names,'x'] = ident
    new_ident = as.character(new_ident$x)
    seur1 = set.ident(seur1, ident.use=new_ident)
}

pdf(paste(args$out, '.tsne.pdf', sep=''), width=12, height=8)
if(args$label == FALSE){
    tsne.plot(seur1, pt.size=.5)
} else{
    tsne.plot(seur1, pt.size=.5, do.label=T)
}
dev.off()
