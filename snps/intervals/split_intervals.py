import re, sys

fn = sys.argv[1]
n = int(sys.argv[2])

prefix = re.search('(.*).intervals', fn).group(1)
print prefix


tot = 0
for line in open(sys.argv[1]):
    if line.startswith('HUMAN'):
        q = map(int, line.rstrip().split(':')[1].split('-'))
        if q[0] > q[1]:
            quit()
        tot += q[1] - q[0]

a = 0
size = 1.*tot/n
cum = 0
i = 0
for line in open(sys.argv[1]):
    if line.startswith('HUMAN'):
        q = line.rstrip().split(':')[1]
        Q = map(int, q.split('-'))
        s = Q[1] - Q[0]
        if cum == 0 or cum >= size:
            try:
                out.close()
            except:
                pass
            out_fn = '%s.%d.intervals' %(prefix, i)
            out = open(out_fn, 'w')
            a += cum
            cum = 0
            i += 1
        out.write(line)
        cum += s
