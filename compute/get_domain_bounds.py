import sys, os
sys.path.append(os.environ['raco'])

def get_domain_bounds(dirname):
    try:
        rmin, rmax = get_parameter(dirname, 'rmin'),\
                get_parameter(dirname, 'rmax')
        nr = get_parameter(dirname, 'n_r')
        domain_bounds = (rmin, rmax)
        ncheby = (nr,)
    except:
        domain_bounds = tuple(get_parameter(dirname, 'domain_bounds'))
        ncheby = tuple(get_parameter(dirname, 'ncheby'))
    return ncheby, domain_bounds
