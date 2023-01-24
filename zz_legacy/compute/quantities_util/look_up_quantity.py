"""
Module to search the lookup table for Rayleigh outputs.

To extract the proper data requires three steps:

     1) read the datafile, e.g. a Shell_Avgs file

           shell_avg = ShellAverage(file, path='')

     2) search the lookup table for the proper indices

           2a) get the lookup table index

                 vr_num = general_lookup(q='v_r')

           2b) get quantity index using the data's lookup table

                 vr_index = shell_avg.lut[vr_num]

     3) extract the data

           vr = shell_avg.data[:,0,vr_index]

There are 3 main routines for searching the table:

     1) general_lookup(q=None, index=None, list=False, alt=False)
     2) shortcut_lookup(quantity, list=False, alt=False, quiet_errors=False)
     3) quantity_in_lut(quantity, alt=False)

     1) Search the table. If you give the index, it returns the variable name.
        If you give the quantity, it returns the index (most often what you want).
        The list keyword will print out the key-value pairs and exit. The alt
        keyword searches the alternate lookup table. For example,
           general_lookup(q='entropy') will return 64
           general_lookup(q='s') will also return 64, see 2) below for why
        while
           general_lookup(index=4, alt=True) will return 'temperature' (which is the entropy)
           general_lookup(index=4) will return 'vp_r' (which is the fluctuating v_r)

     2) Search for the index of a common quantity using a shortcut name. For example
        searching for 'ke' will actually search the table for 'kinetic_energy' and
        's' will search for 'entropy'. This saves on typing and not needing to remember
        the exact variable name in the lookup table. The list keyword will print the
        available shortcut quantities and the alt keyword is the same as before. The
        quiet_errors is more of an internal thing, it suppresses the error output when
        a shortcut is not recognized. For example,
           shortcut_lookup('s') will return 64
           shortcut_lookup('s', alt=True) will return 4

     3) Boolean function that answers "is the given quantity in the
        lookup table?" where quantity is a string like 'v_r' or 'kinetic_energy'.
        It also accepts shortcut names such as 'vr' and 'ke'. For example,
           quantity_in_lut('radial_mke') will return True, mean radial KE is in lookup table
           quantity_in_lut('mrke') will return True since 'mrke' is the shortcut for 'radial_mke'
           quantity_in_lut('mrke', alt=True) will return False, radial_mke is not in older version

R. Orvedahl 7-27-2016

Usage:
    look_up_quantity.py [options]

Options:
    --alt            Use the alternate look up table for "backwards compatability" [default: False]
    --shortcut=<s>   Return the index for shortcut <s> [default: s]
    --quantity=<q>   Return the index for quanityt <q> [default: v_theta]
    --index=<i>      Return quantity with index <i>    [default: 3]
    --list           Print the look up table [default: False]
    --latex=<o>      Write look up table to <o> in LaTeX table format [default:]
"""

from __future__ import print_function
import sys
from collections import OrderedDict
from lut import lut,alt_lut

def quantity_in_lut(quantity, alt=False):
    """
    boolean function: is quantity name in look up table
    """
    q_in_lut = False
    if (alt):
        table = alt_lut # use alternate look up table
    else:
        table = lut

    quantity = quantity.lower()

    # first check the shortcuts
    ind = shortcut_lookup(quantity, alt=alt, quiet_errors=True)
    if (ind is not None):
        q_in_lut = True

    if (not q_in_lut):
        # keep looking and check the table
        for k in table.keys():
            if (quantity == table[k].lower()):
                q_in_lut = True
                break

    return q_in_lut

def general_lookup(q=None, index=None, list=False, alt=False):
    """
    search both shortcuts and the full look up table
    """
    if (q is None and index is None):
        print("\n---ERROR: must specify quantity or index to search the lookup table\n")
        return None

    # only index was given
    if (q is None):
        ind = lookup_full_table(index=index, list=list, alt=alt)
        return ind

    # quantity was given, search shortcuts first
    ind = shortcut_lookup(q, list=list, alt=alt, quiet_errors=True)
    if (ind is not None):
        return ind

    # not found in shortcuts, search full table
    ind = lookup_full_table(q=q, list=list, alt=alt)
    return ind

def lookup_full_table(q=None, index=None, list=False, alt=False):
    """
    given a full quantity name, return the corresponding index or
    given an index, return the corresponding full quantity name
    """
    if (q is None and index is None):
        print("\n---ERROR: must specify quantity or index to search the lookup table\n")
        return None

    if (alt):
        table = alt_lut # use alternate look up table
    else:
        table = lut

    if (list):
        print("\n\tList of index-quantity pairs")
        print("--------------------------------------")
        for k in table.keys():
            print("  {:>26s}\t{:d}".format(table[k],k))
        return None

    # index was given, go get the quantity
    if (index is not None):
        if (index in table.keys()):
           return table[index]
        else:
           return None

    q = q.lower()
    # quantity was given, look up the index
    for k in table.keys():
        if (table[k].lower() == q):
            return k

    print("\n---ERROR: could not find data in look up table, quantity = {}, index = {}\n".format(q, index))
    return None

def shortcut_lookup(quantity, list=False, alt=False, quiet_errors=False):
    """
    use shortcuts to get the index instead of the lengthy exact name
    """
    quantity = quantity.lower()

    shortcuts = OrderedDict()
    shortcuts['ke'] = 'kinetic_energy'
    shortcuts['rke'] = 'radial_ke'
    shortcuts['tke'] = 'theta_ke'
    shortcuts['pke'] = 'phi_ke'
    shortcuts['mrke'] = 'radial_mke'
    shortcuts['mtke'] = 'theta_mke'
    shortcuts['mpke'] = 'phi_mke'

    shortcuts['vr'] = 'v_r'
    shortcuts['vt'] = 'v_theta'
    shortcuts['vp'] = 'v_phi'

    shortcuts['rhovr'] = 'rhov_r'
    shortcuts['rhovt'] = 'rhov_theta'
    shortcuts['rhovp'] = 'rhov_phi'

    shortcuts['s'] = 'entropy'
    shortcuts['dsdr'] = 'entropy_dr'

    shortcuts['cond_flux'] = 'cond_flux_r'
    shortcuts['vol_heat_flux'] = 'vol_heat_flux'
    shortcuts['enth_flux'] = 'enth_flux_radial'
    shortcuts['visc_flux'] = 'visc_flux_r'
    shortcuts['ke_flux'] = 'ke_flux_radial'

    # ---> magnetic quantities <---
    shortcuts['me'] = 'magnetic_energy'

    shortcuts['br'] = 'b_r'
    shortcuts['bt'] = 'b_theta'
    shortcuts['bp'] = 'b_phi'

    # overwrite/add the necessary ones when using the alternate lut
    if (alt):
        shortcuts['dsdr'] = 'gradt_r'  # temperature/entropy are the same
        shortcuts['dtdr'] = 'gradt_r'
        shortcuts['p'] = 'pressure'
        shortcuts['s'] = 'temperature'
        shortcuts['t'] = 'temperature'
        shortcuts['vol_heat_flux'] = 'vol_heating'
        shortcuts['tke'] = 'merid_ke'
        shortcuts['pke'] = 'zonal_ke'
        del shortcuts['mrke'] # there is no equivalent of these
        del shortcuts['mtke'] # in the older version
        del shortcuts['mpke']

    if (list):
        print("\nAvailable shortcuts are:\n")
        for key in shortcuts.keys():
            print("\t{}".format(key))
        print()
        return None

    if (quantity not in shortcuts.keys()):
        if (not quiet_errors):
            print("\n---ERROR: shortcut quantity = {} not recognized\n".format(quantity))
        return None

    q = shortcuts[quantity]
    ind = lookup_full_table(q=q, alt=alt)
    return ind

def write_latex_table(filename, alt=False):
    """
    write look up table to a file in the LaTeX table format
    """
    print("\nwriting latex table...")

    terminate_line = "\n"

    sep = " & "
    end_line = " \\\\" # want '\\' so need to escape it twice

    if (alt):
        table = alt_lut # use alternate look up table
    else:
        table = lut

    mf = open(filename, "w")
    mf.write("Index"+sep+"Quantity Name"+end_line+terminate_line)
    for key in table.keys():
        line = str(key)+sep+str(table[key])+end_line+terminate_line
        mf.write(line)
    mf.close()

    print("saved table to: {}".format(filename))

    return

if __name__ == "__main__":

    from docopt import docopt
    args = docopt(__doc__)

    alt = args['--alt']
    short = args['--shortcut']
    q = args['--quantity']
    ind = int(args['--index'])
    list = args['--list']
    latex = args['--latex']

    if (list):
        lookup_full_table(list=True, alt=alt)
        print()
        ind = shortcut_lookup('', list=True, alt=alt)
        print()
    else:
        print()
        print("shortcut: entered {}, returned index = {}".format(short, shortcut_lookup(short, alt=alt)))
        print("quantity: entered {}, returned index = {}".format(q, lookup_full_table(q=q, alt=alt)))
        print("index   : entered {}, returned quantity = {}".format(ind, lookup_full_table(index=ind, alt=alt)))
        print()

    if (latex != None):
        write_latex_table(latex, alt=alt)

