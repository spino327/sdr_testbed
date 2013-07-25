'''
Copyright (c) 2011, Universidad Industrial de Santander, Colombia
University of Delaware
All rights reserved.

@author: Sergio Pino
@author: Henry Arguello
Website: http://www.eecis.udel.edu/
emails  : sergiop@udel.edu - henarfu@udel.edu
Date   : Dec, 2012
'''

import os

def ensure_dir(f):
    '''
    Ensure that the directory f had been created
    based on code form http://stackoverflow.com/questions/273192/python-best-way-to-create-directory-if-it-doesnt-exist-for-file-write
    '''
    d = os.path.dirname(f)
    if not os.path.exists(d):
        print("i: creating folder " + f + "...")
        os.makedirs(d)