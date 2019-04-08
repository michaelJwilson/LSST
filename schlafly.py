##  https://arxiv.org/pdf/1012.4804.pdf
##  Table 6. 

def get_SchlaflyExt(EBV, filter, Rv=None): 
    ##  7000K source spectrum.
    if Rv == None:
        coeffs = {'u': 4.145, 'g': 3.739, 'r': 2.273, 'i': 1.684, 'z': 1.323, 'y': 1.088} 

        return EBV * coeffs[filter]
        
    else:
        raise UserWarning('Only Rv=3.1 is available.')
