import sys

def get_index(lut, quantity):
    found_quantity = False
    for key in lut.keys():
        if (lut[key].lower() == quantity):
            found_quantity = True
            index = key
    if found_quantity:
        return index
    else:
        print ('The quantity ', quantity, ' does not exist in Rayleigh')
        print ('Exiting ...')
        sys.exit()
