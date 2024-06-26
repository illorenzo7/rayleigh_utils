import subprocess

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
    'Sandybridge',\
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
    'san',\
    'sky_ele',\
    'sky_gpu']

headings = \
    ['Cores',\
    'Total',\
    'Down',\
    'Reserved',\
    'Used',\
    'Free',\
    'Running jobs',\
    'Queued jobs want']

headings_alt = \
    ['cores per node',\
    'nodes total',\
    'nodes down',\
    'nodes reserved',\
    'nodes used',\
    'nodes free',\
    'running jobs',\
    'nodes requested']

node_info = dict({})
for key in keys_alt:
    di_loc = dict({})
    ikey = st.find(key)
    st2 = st[ikey:]
    ihead = 0
    for heading in headings:
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
        head_alt_loc = headings_alt[ihead]
        ihead += 1
        di_loc[head_alt_loc] = num
    node_info[key] = di_loc

for key in keys_alt:
    print(key, node_info[key])
req = dict({})
st2 = st[st.find('requesting'):st.find('using')].split()

count = 0
for subst in st2:
    if subst[:-1] in keys: # remove trailing comma
        req[subst[:-1]] = int(st2[count - 1])
    elif subst == 'Electra':
        nextone = st2[count + 1]
        if nextone == '(B),':
            req['Broadwell (Electra)'] = int(st2[count - 1])
        elif nextone == '(S),':
            req['Skylake'] = int(st2[count - 1])
    elif subst == 'Aitken':
        nextone = st2[count + 1]
        if nextone == '(C),':
            req['Cascadelake'] = int(st2[count - 1])
        elif nextone == '(R)':
            req['ROME'] = int(st2[count - 1])
    count +=1

# get the "using" numbers as well
use = dict({})
st2 = st[st.find('using'):].split()

count = 0
for subst in st2:
    if subst[:-1] in keys: # remove trailing comma
        use[subst[:-1]] = int(st2[count - 1])
    elif subst == 'Electra':
        nextone = st2[count + 1]
        if nextone == '(B),':
            use['Broadwell (Electra)'] = int(st2[count - 1])
        elif nextone == '(S),':
            use['Skylake'] = int(st2[count - 1])
    elif subst == 'Aitken':
        nextone = st2[count + 1]
        if nextone == '(C),':
            use['Cascadelake'] = int(st2[count - 1])
        elif nextone == '(R)':
            use['ROME'] = int(st2[count - 1])
    count +=1

# now print the ratios in a table
def fill_str(stri, lent, char):
    len_loc = len(stri)
    nfill = lent - len_loc
    return stri + char*nfill

char = ' '

headers_care = ['node type', 'requesting', 'using', 'free', 'req. ratio', 'use ratio', 'cores free']
col_buffer = 5
whole_header = ''
column_widths = []
for header in headers:
    if header == headers[-1]:
        column_width = len(header)
    elif header == headers[0]:
        column_width = 25  # needs to be longer due to (Broadwell (Electra))
    else:
        column_width = len(header) + col_buffer
    whole_header += fill_str(header, column_width, char)
    column_widths.append(column_width)

print (whole_header)

# combine both the broadwell types, since they are interchanged freely
# in #PBS requests
rmkey = 'Broadwell (Electra)'
keys.remove(rmkey)
req['Broadwell'] += req[rmkey]
del req[rmkey]
use['Broadwell'] += use[rmkey]
del use[rmkey]
total['Broadwell'] += total[rmkey]
del total[rmkey]

# set number of cores per node
di_ncores = dict({'SandyBridge': 16, 'IvyBridge': 20, 'Haswell': 24, 'Broadwell': 28, 'Skylake': 40, 'Cascadelake': 40, 'ROME':128})

for key in keys:
    nfree = total[key] - use[key]
    ncores = di_ncores[key]
    rowdata = [key, '%i' %req[key], '%i' %use[key], '%i' %(nfree), '%.3f' %(req[key]/total[key]), '%.3f' %(use[key]/total[key]), '%i' %(nfree*ncores)]
    row = ''
    count = 0
    for datum in rowdata:
        row += fill_str(datum, column_widths[count], char)
        count +=1
    print (row)
