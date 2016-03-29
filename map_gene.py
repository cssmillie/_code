import argparse

# Get input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--gene', help='gene id')
parser.add_argument('--genes', help='list of gene ids')
parser.add_argument('--db', help='gene info database', default='/home/unix/csmillie/aviv/db/human.gene_info.txt')
parser.add_argument('--target', help='target list', default='/home/unix/csmillie/aviv/db/human_genes.txt')
args = parser.parse_args()

# Get target genes
target = {}
for line in open(args.target):
    target[line.rstrip()] = 1

# Get query genes
query = {}
if args.gene:
    query[args.gene] = []
if args.genes:
    for line in open(args.genes):
        query[line.rstrip()] = []

# Read gene info
for line in open(args.db):

    if line.startswith('#'):
        continue

    line = line.rstrip().split('\t')
    genes = [line[2]] + line[4].split('|')
    
    qflag = 0
    qgene = ''
    for qgene in genes:
        if qgene in query:
            qflag = 1
            break
    
    if qflag == 0:
        continue
    
    for tgene in genes:
        if tgene in target:
            query[qgene].append(tgene)
    
# Print map
for qgene in query:
    print '%s\t%s' %(qgene, ','.join(query[qgene]))
