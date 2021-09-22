import subprocess

# first get the output of node_stats.sh as a string
st=subprocess.run('node_stats.sh',stdout=subprocess.PIPE)
st=st.stdout.decode("utf-8")

i0 = st.find('hardware type')

keys =\
    ['SandyBridge',\
    'IvyBridge',\
    'Haswell',\
    'Broadwell',\
    'Broadwell (Electra)',\
    'Skylake',\
    'Cascadelake',\
    'ROME']

total = dict({})
for key in keys:
    st2 = st[i0:]
    ikey = st2.find(key)
    st3 = st2[ikey:]
    itot = st3.find('Total')
    st4 = st3[itot:].split()
    total[key] = int(st4[1][:-1]) # get rid of the trailing comma

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

headers = ['node type', 'requesting', 'using', 'free', 'req. ratio', 'use ratio', 'cores free']
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
