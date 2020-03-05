import sys, os

args = sys.argv
nargs = len(args)
fname = None
dirname = None
n = 1000
for i in range(nargs):
    arg = args[i]
    if arg == '-fname':
        fname = args[i+1]
    elif arg == '-dirname':
        dirname = args[i+1]
    elif arg == '-n':
        n = int(args[i+1])

if dirname is None:
    dirname = os.getcwd()

if fname is None:
    files = []
    for name in os.listdir(dirname):
        if 'logfile' in name:
            files.append(name)
    fname = files[-1]

def tail(fname, n):
    f = open(fname, "r")
    """Returns the last `window` lines of file `f` as a list.
    """
    window = n
    if window == 0:
        return []

    BUFSIZ = 1024
    f.seek(0, 2)
    remaining_bytes = f.tell()
    size = window + 1
    block = -1
    data = []

    while size > 0 and remaining_bytes > 0:
        if remaining_bytes - BUFSIZ > 0:
            # Seek back one whole BUFSIZ
            f.seek(block * BUFSIZ, 2)
            # read BUFFER
            bunch = f.read(BUFSIZ)
        else:
            # file too small, start from beginning
            f.seek(0, 0)
            # only read what was not read
            bunch = f.read(remaining_bytes)

        bunch = bunch.decode('utf-8')
        data.insert(0, bunch)
        size -= bunch.count('\n')
        remaining_bytes -= BUFSIZ
        block -= 1
    f.close()
    return ''.join(data).splitlines()[-window:]

print ("dirname: ", dirname)
print("fname: ", fname)

lines = tail(fname, n)
delta_t = []
iters_per_sec = []

for line in lines:
    if "Iteration:" in line:
        split = line.split()
        for i in range(len(split)):
            if split[i] == 'DeltaT:':
                delta_t.append(float(split[i+1]))
            elif split[i] == 'Iter/sec:':
                iters_per_sec.append(float(split[i+1]))

n_final = len(iters_per_sec)
av_dt = np.mean(delta_t)
av_rate = np.mean(iters_per_sec)

print("over last %i iterations" %n_final)
print("avg. time step: %.2f sec" %av_dt)
print("avg. compute rate: %.2f iters/sec" %av_rate)
