import sys
from common import *
from cla_util import *

clas = read_clas(sys.argv)

for key in clas.keys():
    cla = clas[key]
    print (fill_str(key, 25, ' '), cla)

print ("=====================================")
for key in clas_default.keys():
    cla = clas_default[key]
    print (fill_str('default ' + key, 25, ' '), cla)
