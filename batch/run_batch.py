import argparse, ssub

Submitter = ssub.Submitter()
Submitter.p = True

# Input arguments
parser = argparse.ArgumentParser()
parser.add_argument('--seur', help='Seurat file')
parser.add_argument('--field', help='Batch field')
parser.add_argument('--out', help='Output prefix')
args = parser.parse_args()

# Predict cell types
cmd = 'Rscript /home/unix/csmillie/aviv/code/cell_types/predict_types.r --seur %s --labels /home/unix/csmillie/aviv/code/batch/dmap.batch_labels.txt --out %s' %(args.seur, args.out)
Submitter.m = 8
Submitter.submit_and_wait(cmd)

# Get output file
types_fn = '%s.cell_types.txt' %(args.out)

# Correct batch
cmd = 'Rscript /home/unix/csmillie/aviv/code/batch/fix_batch.r --seur %s --field %s --types %s --out %s' %(args.seur, args.field, types_fn, args.out)
Submitter.m = 16
Submitter.q = 'long'
Submitter.submit_and_wait(cmd)

# Get output file
dge_fn = '%s.combat.seur.rds' %(args.out)

# Run seurat
cmd = 'Rscript /home/unix/csmillie/aviv/code/run_seurat.r --seur %s --ident %s' %(args.seur, args.field)
Submitter.submit_and_wait(cmd)

