import subprocess, sys, os
sys.path.append(os.environ['raco'])
from common import fill_str

# first get the output of node_stats.sh as a string
st=subprocess.run('node_stats.sh',stdout=subprocess.PIPE)
st=st.stdout.decode("utf-8")

i1 = st.find('hardware type')
i2 = st.find('GPUs used')
st = st[i1:i2]

keys =\
   ['Broadwell in Pleiades',\
    'Broadwell in Electra',\
    'Cascadelake in Aitken',\
    'Cascadelake in Pleiades',\
    'Haswell',\
    'Ivybridge',\
    'Milan in Cabeus',\
    'Milan in Aitken',\
    'Rome in Aitken',\
    'Rome in Pleiades',\
    #'Sandybridge',\
    'Skylake in Electra',\
    'Skylake in Pleiades']

keys_alt =\
   ['bro',\
    'bro_ele',\
    'cas_ait',\
    'cas_gpu',\
    'has',\
    'ivy',\
    'mil_a100',\
    'mil_ait',\
    'rom_ait',\
    'rom_gpu',\
    #'san',\
    'sky_ele',\
    'sky_gpu']

headers = \
    ['Cores',\
    'Total',\
    'Down',\
    'Reserved',\
    'Used',\
    'Free',\
    'Running jobs',\
    'Queued jobs want']

headers_alt = \
    ['cores per node',\
    'nodes tot.',\
    'nodes down',\
    'nodes res.',\
    'nodes used',\
    'nodes free',\
    'running jobs',\
    'nodes req.']

node_info = dict({})
for key in keys_alt:
    di_loc = dict({})
    ikey = st.find(key)
    st2 = st[ikey:]
    ihead = 0
    print("key=", key)
    for heading in headers:
        i0 = st2.find(heading)
        i0 += len(heading)
        st3 = st2[i0:].split()
        if ihead == 0:
            i0 -= len(heading) # get back to number
            st3 = st2[:i0].split() # get string BEFORE "Cores"
            num = int(st3[-1]) # convert to int
        elif ihead < 7:
            num = int(st3[1][:-1])
            # get rid of the trailing comma and convert to int
        else:
            num = int(st3[0])
        head_alt_loc = headers_alt[ihead]
        ihead += 1
        di_loc[head_alt_loc] = num
    node_info[key] = di_loc

#for key in keys_alt:
#    print(key, node_info[key])

# combine the broadwells 
node_info['bro_tot'] = dict({})
keys.insert(2, 'Broadwell (Total)')
keys_alt.insert(2, 'bro_tot')

for key in node_info['bro'].keys():
    node_info['bro_tot'][key] = \
            node_info['bro'][key] + node_info['bro_ele'][key]

# now print the ratios in a table
char = ' '

headers_care = \
    ['node type',
    'req. ratio',\
    'use ratio',\
    'down ratio',\
    'nodes tot.',\
    'nodes down',\
    'nodes res.',\
    'nodes used',\
    'nodes free',\
    'nodes req.']

col_buffer = 5
whole_header = ''
column_widths = []
for header in headers_care:
    if header == headers[-1]:
        column_width = len(header)
    elif header == headers[0]:
        column_width = 25  # needs to be longer due to (Broadwell (Electra))
    else:
        column_width = len(header) + col_buffer
    whole_header += fill_str(header, column_width, char)
    column_widths.append(column_width)

print (whole_header)

for key in keys_alt:
    di_loc = node_info[key]
    ncores = di_loc['cores per node']
    ntot = di_loc['nodes tot.']
    nuse = di_loc['nodes used']
    ndown = di_loc['nodes down']
    nreq = di_loc['nodes req.']

    rowdata = []
    for header in headers_care:
        if header == 'node type':
            rowdata += [key]
        elif header == 'req. ratio':
            rowdata += ['%.3f' %(nreq/ntot)]
        elif header == 'use ratio':
            rowdata += ['%.3f' %(nuse/ntot)]
        elif header == 'down ratio':
            rowdata += ['%.3f' %(ndown/ntot)]
        else:
            rowdata += ['%i' %di_loc[header]]
    #rowdata = [key, '%i' %req[key], '%i' %use[key], '%i' %(nfree), '%.3f' %(req[key]/total[key]), '%.3f' %(use[key]/total[key]), '%i' %(nfree*ncores)]
    row = ''
    count = 0
    for datum in rowdata:
        row += fill_str(datum, column_widths[count], char)
        count +=1
    print (row)
