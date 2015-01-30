# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 23:45:12 2015

@author: Danielbraga
"""
#FUNCTION: 
#def function():
handle = open("uniprot.txt").readlines()
function='CC   '
uniprotid='AC' 
listfunction=[]
uniprot=[]
lista=[]

for i in handle:
    if function in i:
        y=i.split(" ")             
        listfunction.append(i)
#
#for j in handle:
#    if uniprotid in j:
#        uniprot.append(j)
#x=len(listname)        
#for k in range(x):
#    lista.append(((uniprot[k]+" "+listname[k])))
#return lista
##print(review())
