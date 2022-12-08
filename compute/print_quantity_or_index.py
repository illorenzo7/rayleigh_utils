"""
Author: Loren Matilsky
Created: 09/13/21
Wrapper for Ryan's parse_quantity

This script prints the Rayleigh index or variable name associated with a
given quantity.

Takes as an argument a user-specified diagnostic quantity or 
quantity index, and translates between the two."""

import sys, lut
from common import make_bold
quantity_or_index = sys.argv[1]

index, quantity = lut.parse_quantity(quantity_or_index)
if index is None or quantity is None:
    print ('The index or quantity', make_bold(quantity_or_index), 'does not exist in Rayleigh!')
else:
    print ('The index', make_bold(str(index)), 'corresponds to the quantity', make_bold(quantity))
