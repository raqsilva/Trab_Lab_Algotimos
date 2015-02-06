# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 14:21:20 2015

@author: Danielbraga
"""
from Bio import SeqIO 
from Bio import Entrez
from Bio.Blast import NCBIWWW,NCBIXML

import shutil#moving files
import os.path#cheking files in path
import urllib
#from Uniprot_Parser import * #parsing uniprot text file
from Bio.SeqIO import UniprotIO #parsing uniprot xml file
import pandas
import numpy as np
import urllib.request
import os
def allhits():
    lista=[]
    handle = open("../res/blast_without_note/match/matches.txt").readlines()    
    first = 'sp|'
    last = '|'
    gi=[]
    xml='['
    
    
        
    for q in range(len(handle)):
        lista.append([])
        x=handle[q].split()
        for i in range (len(x)):
            m=x[i]
            if 'sp|' in m:    
                try:
                    start = m.rindex( first ) + len( first )
                    end = m.rindex( last, start )
                    novo= m[start:end]
                    lista[q].append(novo)
                except:
                    pass
          
    
    for f in handle:
        n = f[f.index(xml) + len(xml):] 
        g = n.split('.xml', 1)[0]  
        gi.append(g[1:])
    
    for p in range(len(lista)):
        lista[p].append(gi[p])
    return lista
file = open("../res/blast_without_note/match/function/allhits.txt",'w')
lista=allhits()         
for i in range(len(lista)):
        file.write("%s\n" % lista[i])
file.close()