'''
Copyright (c) 2012, Universidad Industrial de Santander, Colombia
University of Delaware
All rights reserved.

@author: Sergio Pino
@author: Henry Arguello
Website: http://www.eecis.udel.edu/
emails  : sergiop@udel.edu - henarfu@udel.edu
Date   : Dec, 2012
'''

def readProperties(cfile):

    properties = dict()
    
    try:
        confFile = open(cfile, "r")
    
        for line in confFile:
            line = line.strip()
    
            if len(line) > 0:
                info = line.strip().split(":", 1)
    
                if len(info) != 2:
                    raise MemoryError("e: Unable to retrieve conf information")
                
                properties[info[0]] = info[1]
    
    except Exception, e:
        print(e)
        confFile.close()
    
    return properties
