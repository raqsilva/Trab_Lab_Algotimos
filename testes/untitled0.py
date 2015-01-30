# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 23:45:12 2015

@author: Danielbraga
"""
#name , status of review and size
def review():
    handle = open("uniprot.txt").readlines()
    status='ID   '
    aa='AA.'
    res=[]
    
    
    for i in handle:
        if status and aa in i:  
            res.append(i) 
    return res
    
print(review())
