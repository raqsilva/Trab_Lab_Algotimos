# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 23:45:12 2015

@author: Danielbraga
"""


handle = open("uniprot.txt").readlines()
uniprotid='AC   '
mystr = 'SUBCELLULAR LOCATION:'
lista=[]
l=[]
for i in handle:
    if mystr in i:
        y=i.split(" ")            
        lista.append(y[6] )       

for j in handle:
    if uniprotid in j:
        l.append(j)

