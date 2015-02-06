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

handle = open("../res/blast_without_note/match/matches.txt").readlines()

gimatch=[]
gi=[]
hit=[]
s2='sp|'
xml='['
for i in handle:
    n = i[i.index(xml) + len(xml):] 
    g = n.split('.xml', 1)[0]  
    gi.append(g)
for J in handle:
    s3 = J[J.index(s2) + len(s2):] 
    sP = s3.split('|', 1)[0] 
    sp=sP[:6]    
    hit.append(sp)
for k  in range(len(gi)):
    gimatch.append(gi[k]+' '+hit[k])    
        
        
for w in range(len(hit)):  
    site = urllib.request.urlopen("http://www.uniprot.org/uniprot/"+hit[w]+".txt")
    data = site.readlines()
    nome=hit[w]
    file = open("../res/blast_without_note/match/function/"+nome+'.txt',"wb") #open file in binary mode
    file.writelines(data)
    file.close()