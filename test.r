library(optparse)
if(interactive()){
    args = list()
    args$map = '~/aviv/summary/test.map'
    args$ming = 1000
    args$npts = 10
    args$pseudocount = 0
    args$out = 'test'
} else{
option_list = list(make_option('--map', help='input map (1=id, 2=folder, 3=pattern)'),
                   make_option('--ming', help='genes per cell cutoff', type='integer'),
                   make_option('--npts', help='number of points to calculate for each plot'),
                   make_option('--pseudocount', default=0, type='integer'),
                   make_option('--out', help='output prefix')
                   )
args = parse_args(OptionParser(option_list=option_list))
}

multiplot = function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    require(grid)
    plots = c(list(...), plotlist)
    numPlots = length(plots)
    if (is.null(layout)) {
        layout = matrix(seq(1, cols * ceiling(numPlots/cols)),
                        ncol = cols, nrow = ceiling(numPlots/cols))
    }
    if (numPlots==1) {
        print(plots[[1]])
    } else {
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        for (i in 1:numPlots) {
            matchidx = as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}


data = read.table(args$map, stringsAsFactors=F)

g = c()
x = c()
y.h = c()
y.m = c()
y.t = c()

for(i in 1:nrow(data)){

    id = data[i,1]
    folder = data[i,2]
    pattern = data[i,3]

    h.fn = paste('~/aviv/data/', folder, '/dge/human.dge.txt.gz', sep='')
    m.fn = paste('~/aviv/data/', folder, '/dge/mouse.dge.txt.gz', sep='')

    h.dge = read.table(paste('~/aviv/data/', folder, '/dge/human.dge.txt.gz', sep=''), header=T, row.names=1, sep='\t')
    h.use = grep(pattern, colnames(h.dge))
    h.dge = as.matrix(h.dge[,h.use])
    
    m.dge = read.table(paste('~/aviv/data/', folder, '/dge/mouse.dge.txt.gz', sep=''), header=T, row.names=1, sep='\t')
    m.use = grep(pattern, colnames(m.dge))
    m.dge = as.matrix(m.dge[,m.use])

    total_reads = 0
    if(ncol(h.dge) > 1){total_reads = total_reads + sum(h.dge + args$pseudocount)}
    if(ncol(m.dge) > 1){total_reads = total_reads + sum(m.dge + args$pseudocount)}

    max_reads = log10(total_reads)
    min_reads = log10(10**5)
    x.i = seq(min_reads, max_reads, length.out=args$npts)
    g = c(g, rep(id, length(x.i)))

    ncoll = 0
    
    if(ncol(h.dge) > 1){
        h.coll = length(unique(sapply(strsplit(colnames(h.dge), '\\.'), '[', 1)))
        ncoll = h.coll
        h.prob = as.vector(h.dge + args$pseudocount)/total_reads
        for(i in x.i){
            num_reads = 10**i
            h.tmp = matrix(rmultinom(1, num_reads, h.prob), nrow=nrow(h.dge), ncol=ncol(h.dge))
            h.cells = sum(colSums(h.tmp > 0) >= args$ming)
            print(c('human', h.cells))
            y.h = c(y.h, h.cells/h.coll)
        }
    } else {
        y.h = c(y.h, rep(0, length(x.i)))
    }

    if(ncol(m.dge) > 1){
        m.coll = length(unique(sapply(strsplit(colnames(m.dge), '\\.'), '[', 1)))
        ncoll = m.coll
        m.prob = as.vector(m.dge + args$pseudocount)/total_reads
        for(i in x.i){
            num_reads = 10**i
            m.tmp = matrix(rmultinom(1, num_reads, m.prob), nrow=nrow(m.dge), ncol=ncol(m.dge))
            m.cells = sum(colSums(m.tmp > 0) >= args$ming)
            print(c('mouse', m.cells))
            y.m = c(y.m, m.cells/m.coll)
        }
    } else {
        y.m = c(y.m, rep(0, length(x.i)))
    }
    x = c(x, x.i/ncoll)
}
y.t = y.h + y.m


library(ggplot2)
