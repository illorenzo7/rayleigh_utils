import sys
from common import *

args = sys.argv[2:]
clas = read_clas(args)

for key in clas.keys():
    cla = clas[key]
    print (fill_str(key, 25, ' '), cla.val)

print ("=====================================")
for key in clas_default.keys():
    cla = clas_default[key]
    print (fill_str(key, 25, ' '), cla.val)
