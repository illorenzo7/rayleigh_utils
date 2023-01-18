# Author: Loren Matilsky
# Created: well before 05/06/2019
# a string of text (either array, list, int, float etc.) into a the
# appropriate Python data type/value
import numpy as np

def how_to_treat_numberstring(string):
    int_chars = ['+', '-', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    float_chars = int_chars + ['.', 'e', 'd', 'E', 'D']
    bool_strings = ['true', 'false', 't', 'f', '.true.', '.false.']
    if string.lower() in bool_strings:
        return 'bool'
    elif ',' in string:
        return 'array'
    elif all (char in int_chars for char in string):
        return 'int'
    elif all (char in float_chars for char in string):
        if not string[0] in (int_chars + ['.']): # this must be a string, like d10.0
            return 'string'
        elif string[-2:] == '.0': # convert things like 1e3 --> '1000.0'
            # to ints
            return 'int'
        else:
            return 'float'    
    else: # just a string!
        return 'string'

def trim_arraystring(string):
    # Ignore possible enclosing braces. Note that this will screw up things
    # we ever want to deal with multidimensional arrays
    trimmed_string = string.replace('[', '')
    trimmed_string = trimmed_string.replace('(', '')
    trimmed_string = trimmed_string.replace('{', '')
    trimmed_string = trimmed_string.replace(']', '')
    trimmed_string = trimmed_string.replace(')', '')
    trimmed_string = trimmed_string.replace('}', '')    
    return trimmed_string

def string_to_bool(string):
    true_strings = ['true', 't', '.true.']
    false_strings = ['false', 'f', '.false.']
    if string.lower() in true_strings:
        return True
    elif string.lower() in false_strings:
        return False
    
def string_to_number_or_array(string):
    number_type = how_to_treat_numberstring(string)
    
    if number_type == 'bool':
        return string_to_bool(string)
        
    elif number_type == 'float':
        string = string.replace('d', 'e')
        string = string.replace('D', 'e')
        return float(string)
    
    elif number_type == 'int':
        return int(string)
    
    elif number_type == 'array':
        mylist = []
        string = trim_arraystring(string)
        string_list = string.split(',')
        for substring in string_list:
            member_number_type = how_to_treat_numberstring(substring)
            if member_number_type == 'bool':
                member_val = string_to_bool(substring)          
            elif member_number_type == 'float':
                substring = substring.replace('d', 'e')
                substring = substring.replace('D', 'e')
                member_val = float(substring)
            elif member_number_type == 'int':
                member_val = int(substring)
            mylist.append(member_val)
        
        return np.array(mylist)
    elif number_type == 'string':
        return string
