"""
Author: Loren Matilsky
Created: 03/27/2017
    This script contains the use_lut(string) which translates a given string     to the corresponding variable index specified in diagnostics_base.mod 

    format:
        index = use_lut(string) 

    e.g., use_lut('entropy') would produce int(64)
"""
from look_up_quantity import *

def use_lut(string):
    return general_lookup(q=string)
