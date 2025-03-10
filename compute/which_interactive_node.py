import subprocess, sys, os
sys.path.append(os.environ['raco'])
from common import fill_str

# see what type of interactive node we're after
searchparameter = 'interactive_quote_unquote'
args = sys.argv[1:]
jobid = None
i = 0
for arg in args:
    if arg == '--stdin':
        searchparameter = 'STDIN'
    elif arg == '--id':
        jobid = args[i+1]
    i += 1

# first get the output of node_stats.sh as a string
st=subprocess.check_output(['qstat', '-u', 'lmatilsk'])
st=st.decode('utf-8')

lines = st.split('\n')
the_line = None
keep_searching = True
for line in lines:
    if keep_searching:
        if searchparameter in line:
            the_line = line
    if not jobid is None:
        if jobid in line:
            keep_searching = False

if the_line is None:
    print ("No interactive session is currently being run.")
    print ("To run one, type ")
    if searchparameter == 'interactive_quote_unquote':
        print ("interactivenode [node_type] [number_of_nodes_desired]")
        print ("e.g.")
        print ("interactivenode has 16")
    elif searchparameter == 'STDIN':
        print ("devel [node_type] [cores]")
        print ("e.g.")
        print ("devel has 256")
else:
    where = the_line.find('.pbspl1')
    job_number = the_line[:where]
    st2=subprocess.check_output(['qstat', '-f', job_number])
    st2=st2.decode('utf-8')
    running = False
    if 'exec_host' in st2:
        running = True

    if not running:
        print ("Interactive job", job_number, "has been requested but isn't running yet. Please wait.")
    else:
        where = st2.find("exec_host = ")
        st3 = st2[where+12:]
        where = st3.find("/0")

        the_node = st3[:where]
        print("the exec node for interactive job", job_number, "is")
        print(the_node)
