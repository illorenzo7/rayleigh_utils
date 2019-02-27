import sys
from string_to_num import string_to_number_or_array

def get_parameter(dirname, parameter):
    f = open(dirname + '/main_input')
    lines = f.readlines()
    n = len(lines)
    try:
        for i in range(n):
            if (parameter in lines[i] and '=' in lines[i]):
                line = lines[i]
        line = line[:] # test if line was assigned
    except:
        print('Note: ' + parameter + ' was not')
        print('specified in run: ' + dirname + '. ')
        print('Assigning the arbitrary value 100')
        return(100)
    
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
