import argparse

# get command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('--map', help='mapping file')
parser.add_argument('--dge', default='', help='input dge')
parser.add_argument('--regex', help='use cells that match regex', default='')
parser.add_argument('--ident', help='identity field separator', default='1')
parser.add_argument('--vcut', help='variable genes cutoff (diffCV)', type=float, default=-1)
parser.add_argument('--step', choices=['all', 'merge', 'tsne', 'cluster', 'markers', 'plots'], default='all')
parser.add_argument('--human', default=False, action='store_true')
parser.add_argument('--mouse', default=False, action='store_true')
parser.add_argument('--out', help='output prefix')
parser.add_argument('--ming', help='comma-separated list of ming', default='500,750,1000')
parser.add_argument('--minc', help='comma-separated list of minc', default='25,50')
parser.add_argument('--retry', help='number of retry attempts', default=0, type=int)
parser.add_argument('--test', default=False, action='store_true')
args = parser.parse_args()

# fix command line arguments
if args.human == False and args.mouse == False:
    quit('error: must specify --human or --mouse')
mings = map(int, args.ming.split(','))
mincs = map(int, args.minc.split(','))
if args.step in ['all', 'merge', 'tsne'] and args.vcut == -1:
    quit('error: must specify --vcut with --step %s' %(args.step))

import os, re, ssub
import numpy as np

# enter pipeline at any step
def check_name(name):
    steps = 'all merge tsne cluster markers plots'.split()
    if steps.index(args.step) <= steps.index(name):
        return True
    else:
        return False

def mkchdir(dname):
    if not os.path.exists(dname):
        os.mkdir(dname)
    os.chdir(dname)
    return True

def check_outs(cmds, outs, count=0):
    I = []
    if count == 0 and len(outs) == 0:
        return range(len(cmds))
    for i,out in enumerate(outs):
        if not os.path.exists(out) or os.stat(out).st_size == 0:
            I.append(i)
    return I

# initialize job submitter
Submitter = ssub.Submitter()
if args.test:
    Submitter.p = True

# run cmd & check output
def run_cmds(name, cmds, outs, q='short', m=16):
    # check step
    if check_name(name) == False:
        return True
    # initialize submitter
    Submitter.q = q
    Submitter.m = m
    # print tasks
    if Submitter.p == True:
        indices = check_outs(cmds, outs, count=0)
        if len(indices) > 0:
            print '\n'.join([cmds[i] for i in indices])
        return True
    # submit tasks
    count = 0
    while True:
        # check retry attempts
        if (count - 1) > args.retry:
            quit('Submission count exceeded --retry')
        # submit unfinished tasks
        indices = check_outs(cmds, outs, count=count)
        if len(indices) > 0:
            cmds = [cmds[i] for i in indices]
            if len(outs) > 0:
                outs = [outs[i] for i in indices]
            print 'Submitting with memory %d and queue %s:\n' %(Submitter.m, Submitter.q) + '\n'.join(cmds) + '\n\n'
            Submitter.submit_and_wait(cmds)
            # update retry variables
            Submitter.m += 8
            count += 1
        else:
            break
    return True


# step 1: merge dges
# merge dges
tmp = '%s.dge.txt' %(args.out)
cmd1 = 'Rscript ~/aviv/code/merge_dges.r --map %s --ming 500 --out %s' %(args.map, tmp)
run_cmds('merge', [cmd1], [tmp], m=16)
# gzip dge
dge = os.path.abspath('%s.gz' %(tmp))
cmd1 = 'gzip %s' %(tmp)
run_cmds('merge', [cmd1], [dge])

if args.dge != '':
    dge = os.path.abspath(args.dge)


# step 2: run seurat
mkchdir('../tsne')
# run seurat
seurs = []
cmds = []
outs = []
logs = []
for ming in mings:
    for minc in mincs:
        out = '%s.%d.%d' %(args.out, ming, minc)
        seurs += [os.path.abspath('%s.seur.rds' %(out))]
        if args.regex != '':
            cmds += ['Rscript ~/aviv/code/run_seurat.r --dge %s --ming %d --minc %d --regex %s --ident %s --out %s --vcut %f' %(dge, ming, minc, args.regex, args.ident, out, args.vcut)]
        else:
            cmds += ['Rscript ~/aviv/code/run_seurat.r --dge %s --ming %d --minc %d --ident %s --out %s --vcut %f' %(dge, ming, minc, args.ident, out, args.vcut)]
        outs += [os.path.abspath('%s.tsne2.pdf' %(out))]
        logs += [os.path.abspath('%s.log.rds' %(out))]
run_cmds('tsne', cmds, outs, q='long', m=16)


# step 3: cluster tsne
mkchdir('../cluster')
# run graph clustering
cmds = []
clusts = []
for i in range(len(seurs)):
    seur = seurs[i]
    log = logs[i]
    prefix = re.search('(.*).seur.rds', os.path.basename(seur)).group(1)
    for nnk in [10, 25, 50, 100]:
        out = '%s.imap.%d' %(prefix, nnk)
        cmds += ['Rscript ~/aviv/code/cluster_graph.r --seur %s --log %s --pca --nnk %d --type infomap --out %s' %(seur, log, nnk, out)]
        clusts += [os.path.abspath('%s.clust.rds' %(out))]
run_cmds('cluster', cmds, clusts, q='long', m=32)


# step 4: markers
mkchdir('../markers')
# find all markers
cmds = []
plot = []
marks = []
for clust in clusts:
    out = re.search('(.*).clust.rds', os.path.basename(clust)).group(1)
    seur = re.search('(.*).imap', clust).group(1) + '.seur.rds'
    seur = re.sub('/cluster/', '/tsne/', seur)
    cmds += ['Rscript ~/aviv/code/find_markers.r --seur %s --clust %s --out %s' %(seur, clust, out)]
    marks += [os.path.abspath('%s.markers.roc.txt' %(out))]
run_cmds('markers', cmds, marks, q='long', m=16)


# step 5: plots
mkchdir('../genes')
# plot genes & cell types
cmds = []
for clust in clusts:
    out = re.search('(.*).clust.rds', os.path.basename(clust)).group(1)
    seur = re.search('(.*).imap', clust).group(1) + '.seur.rds'
    seur = re.sub('/cluster/', '/tsne/', seur)
    if args.mouse == True:
        cmds += ['Rscript ~/aviv/code/cell_types/predict_types.r --seur %s --mouse --out %s' %(seur, out)]
        cmds += ['Rscript ~/aviv/code/plot_genes.r --seur %s --clust %s --genes ~/aviv/db/markers/essential.m.txt --hmap --out %s.essential' %(seur, clust, out)]
        cmds += ['Rscript ~/aviv/code/plot_matrix.r --seur %s --clust %s --matrix ~/aviv/db/immgen/immgen.mm10.txt --hmap --out %s.immgen' %(seur, clust, out)]
    if args.human == True:
        cmds += ['Rscript ~/aviv/code/cell_types/predict_types.r --seur %s --out %s' %(seur, out)]
        cmds += ['Rscript ~/aviv/code/plot_genes.r --seur %s --clust %s --genes ~/aviv/db/markers/essential.h.txt --hmap --out %s.essential' %(seur, clust, out)]
        cmds += ['Rscript ~/aviv/code/plot_matrix.r --seur %s --clust %s --matrix ~/aviv/db/dmap/output/dmap_expression_mapped.tsv --hmap --out %s.dmap' %(seur, clust, out)]
run_cmds('plots', cmds, [], q='long', m=16)
