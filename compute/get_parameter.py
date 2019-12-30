import sys
from string_to_num import string_to_number_or_array

def get_parameter(dirname, parameter):
    f = open(dirname + '/main_input')
    lines = f.readlines()
    n = len(lines)
    try:
        for i in range(n):
            if (parameter in lines[i] and '=' in lines[i] and \
                    lines[i][0] != '!'):
                line = lines[i]
        line = line[:] # test if line was assigned
    except:
#        print('Note: ' + parameter + ' was not')
#        print('specified in run: ' + dirname + '. ')
#        print('exiting ...')
        if parameter == 'magnetism' or parameter == 'use_extrema':
            return False # if magnetism wasn't specified, it is "False"
        else:
            raise Exception('The parameter ' + parameter + ' was not\n' +\
                            'specified in run: ' + dirname + '. \n' +\
                            'exiting ....\n')
    
    # Make line lowercase
    line = line.lower()
    
    # Remove spaces and newline character (at the end of each line)
    line = line.replace(' ', '')
    line = line.replace('\n', '')

    equals_index = line.index('=') # find where the actual number
        # or array starts (should be after the equals sign)
    num_string = line[equals_index + 1:]
    if '!' in num_string: # there was a comment after the equals statement
                        # throw it away!
        excl_index = num_string.index('!')
        num_string = num_string[:excl_index]
    return (string_to_number_or_array(num_string))
