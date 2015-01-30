# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 23:45:12 2015

@author: Danielbraga
"""

def proteinlocation():
    handle = open("uniprot.txt").readlines()
    uniprotid='AC   '
    
    mystr = 'SUBCELLULAR LOCATION:'
    location=[]
    uniprot=[]
    lista=[]
    
    for i in handle:
        if mystr in i:
            y=i.split(" ")            
            location.append(y[6] )
    
    for j in handle:
        if uniprotid in j:
            uniprot.append(j)
    x=len(location)        
    for k in range(x):
        lista.append(((uniprot[k]+" "+location[k])))
        
    return lista
    
x=(proteinlocation())