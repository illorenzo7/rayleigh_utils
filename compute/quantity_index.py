"""
Author: Loren Matilsky
Created: 03/23/2017
Last Modified: 11/18/2018

This script prints the Rayleigh index or variable name associated with a
given quantity.

Takes as an argument a user-specified diagnostic quantity or 
quantity index, and translates between the two. For example, typing
    python quantity_index entropy   
would produce the output '501' on the command line.
    
Conversely, typing
    python quantity_index 64 -old 
would produce the output 'entropy' on the command line.
"""

import sys, lut
from get_index import get_index

lut_current = lut.lut
        
quantity_or_index = sys.argv[1]

try:
    index = int(quantity_or_index) # For '1' or '2', this will be successful
                                # for 'v_r,' this will fail
    print ('The index ', index, ' corresponds to the quantity ',\
           lut_current[index])
except:  # quantity_or_index was NOT '1,' '2,' ... must be quantity string!
    print ('The quantity ', quantity_or_index, ' has index ',\
           get_index(lut_current, quantity_or_index))
