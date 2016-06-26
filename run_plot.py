import argparse, re

# Get command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--seur', help='seurat file (rds)', required=True)
parser.add_argument('--ident', help='cluster file (rds)', default='')
parser.add_argument('--human', help='use human databases', default=False, action='store_true')
parser.add_argument('--mouse', help='use mouse databases', default=False, action='store_true')
parser.add_argument('--single', help='make single plots', default=False, action='store_true')
args = parser.parse_args()

# Get output prefix
prefix = re.search('(.*).seur.rds', args.seur.split('/')[-1]).group(1)

# Construct genes commands
def genes_cmd(seur, ident, genes, out, combine='sum'):
    cmd = 'Rscript /home/unix/csmillie/aviv/code/plot_genes.r --seur %s --genes %s --tsne --out %s --combine %s' %(seur, genes, out, combine)
    return cmd

# Construct matrix commands
def matrix_cmd(seur, ident, genes, out, single):
    cmd = 'Rscript /home/unix/csmillie/aviv/code/plot_matrix.r --seur %s --matrix %s --out %s' %(seur, genes, out)
    if single:
        cmd += ' --tsne_i'
    else:
        cmd += ' --tsne --tsne_k 25'
    if ident:
        cmd += ' --ident %s' %(ident)
    return cmd

if args.human:

    # DMAP
    out = '%s.dmap' %(prefix)
    genes = '/home/unix/csmillie/aviv/db/dmap/output/dmap_expression_mapped.tsv'
    print matrix_cmd(args.seur, args.ident, genes, out, args.single)
    
if args.mouse:

    # Immgen
    out = '%s.immgen' %(prefix)
    genes = '/home/unix/csmillie/aviv/db/immgen/immgen.mm10.txt'
    print matrix_cmd(args.seur, args.ident, genes, out, args.single)

# Plot single markers

combine = ['essential', 'surface', 'g1s', 'g2m', 'tcells']
singles = ['surface', 'bcells']

for name in combine:
    out = '%s.%s' %(prefix, name)
    if args.human:
        genes = '/home/unix/csmillie/aviv/db/markers/%s.h.txt' %(name)
    if args.mouse:
        genes = '/home/unix/csmillie/aviv/db/markers/%s.m.txt' %(name)
    print genes_cmd(args.seur, args.ident, genes, out, combine='sum')

for name in singles:
    out = '%s.%s' %(prefix, name)
    if args.human:
        genes = '/home/unix/csmillie/aviv/db/markers/%s.h.txt' %(name)
    if args.mouse:
        genes = '/home/unix/csmillie/aviv/db/markers/%s.m.txt' %(name)
    print genes_cmd(args.seur, args.ident, genes, out, combine='none')

