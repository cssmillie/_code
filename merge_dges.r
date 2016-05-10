# This script merges multiple DGEs into a single DGE

# Read input arguments
library(optparse)
option_list = list(make_option('--map', help='input mapping file (1=project, 2=path, 3=pattern)'),
                   make_option('--ming', help='genes per cell cutoff', type='integer'),
                   make_option('--out', help='output file')
                   )
args = parse_args(OptionParser(option_list=option_list))

library(plyr)
library(stringr)

# Extract variables
map = read.table(args$map, stringsAsFactors=F)

# Get list of DGEs
dges = list()
for(i in 1:nrow(map)){
    name = map[i,1]
    path = map[i,2]
    pattern = map[i,3]
    data = read.table(path, header=T, row.names=1, sep='\t')
    dim1 = dim(data)
    data = data[,grep(pattern, colnames(data))]
    dim2 = dim(data)
    print(c(name, dim1, dim2))
    dges[[name]] = data
}

# Relabel colnames to prevent collisions
for(i in 1:length(dges)){colnames(dges[[i]]) = paste(map[i,1], colnames(dges[[i]]), sep='.')}

merge_dges = function(dges){
    # merge a list of dges into a single dge
    rows = sort(unique(unlist(llply(dges, rownames)))) # get rownames across all dges
    cols = sort(unique(unlist(llply(dges, colnames)))) # get colnames across all dges
    dge = data.frame(matrix(0, nrow=length(rows), ncol=length(cols))) # initialize dge matrix
    rownames(dge) = rows
    colnames(dge) = cols
    for(i in 1:length(dges)){dge[rownames(dges[[i]]), colnames(dges[[i]])] = dges[[i]]} # add all dges
    return(dge)
}

dge = merge_dges(dges)

# Filter DGE by gene count
dge = dge[, colSums(dge > 0) >= args$ming]

# Write DGE to outfile
write.table(dge, args$out, quote=F, sep='\t')
