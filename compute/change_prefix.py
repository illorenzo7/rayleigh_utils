import sys, os

rootdir, pref1, pref2 = sys.argv[1:]

for fname in os.listdir(rootdir):
    fname_orig = os.path.join(rootdir, fname)
    fname = fname_orig.replace(pref1, pref2)
    os.rename(fname_orig, fname)
