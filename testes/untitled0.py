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
#  handle = open("uniprot.txt").readlines()
#    status='ID   '
#    aa='AA.'
#    res=[]
#    
#    
#    for i in handle:
#        if status and aa in i:  
#            res.append(i) 

blast=[]    
for file in os.listdir("../res/blast_without_note/match/function/"):
    if file.endswith(".txt"):
           blast.append(file)

lista=[]
match=[]
for j in range(len(blast)):
    nome=blast[j]
    gi = nome.replace(".txt","")
    func='CC   -!- FUNCTION:'
    file = open("../res/blast_without_note/match/function/"+nome,'r') 
    data = file.readlines()
    for i in data:
        if func in i:
            function=i.replace('CC   -!- FUNCTION:', '')
            lista.append('Gi: '+gi+'  '+'Possivel função:  '+function)
#function to se the proteins that have possible note:   
for n in range(len(blast)):
    name=blast[n]
    gis = name.replace(".txt","")
    for k in lista:
        if gis in k:
            match.append(gis)