def get_parameter(dirname,parameter):
    '''
    This function searches the 'main_input' file in a Rayleigh run
    (the run associated with 'dirname') and returns the associated 
    simulation parameter. 

    For example, get_parameter('Gastine','rmin') will return the 
    associated minimum radius of the run 'Gastine'. 
    '''
    import sys
    f = open(dirname + '/main_input')
    lines = f.readlines()
    n = len(lines)
    try:
        for i in range(n):
            if (parameter in lines[i]):
                line = lines[i]
        line = line[:]
    except:
        print('Note: ' + parameter + ' was not')
        print('specified in run: ' + dirname + '. ')
        print('Assigning the arbitrary value 1.0')
        return(1.0)
    equals_index = line.index('=')
    i = equals_index + 1
    while (not is_part_of_number(line[i])):
        i = i + 1
    number_start = i + 0 #avoid copying the reference instead of the value
    j = number_start + 0
    while (is_part_of_number(line[j])):
        j = j + 1
    number_end = j + 0
    string_number = line[number_start:number_end]
    list_number = list(string_number)
    len_number = len(string_number)
    for i in range(len_number):
        if (list_number[i] == 'd' or list_number[i] == 'D'):
            list_number[i] = 'e' # Python floats don't like "D"!
    string_number = ''.join(list_number)
    return float(string_number)

def is_part_of_number(char):
    if (ord('0') <= ord(char) and ord(char) <= ord('9')):
        return True
    elif (char == 'd' or char == 'D' or char == 'e' or char == 'E'):
        return True
    elif (char == '+' or char == '-' or char == '.'):
        return True
    else:
        return False

