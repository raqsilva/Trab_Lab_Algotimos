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

blast=[]    
for file in os.listdir("../res/blast_with_note"):
    if file.endswith(".xml"):
        blast.append(file)
E_VALUE_THRESH = 0.05
lista=[]
for i in range(len(blast)):
    lista.append([])
    lista[i].append(blast[i])
    result_handle = open("../res/blast_with_note/"+blast[i])
    blast_record = NCBIXML.read(result_handle)
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
             if hsp.expect < E_VALUE_THRESH:
                 lista[i].append(alignment.title)
                 lista[i].append(alignment.length)
                 lista[i].append(hsp.expect)
save_file = open('nomatches.txt', "w")
for i in range(len(lista)):                
    if len(lista[i])<2:
        save_file.write(str(lista[i])+'\n')
save_file.close()
#        #moving the file to another directory
path=os.getcwd()
src = path+"/"+'nomatches.txt' #source folder
dst = "../res/blast_with_note/nomatch/"#destination folder
shutil.move(src, dst)
                 
save_file = open('matches.txt', "w")
for i in range(len(lista)):                
    if len(lista[i])>2:
        save_file.write(str(lista[i])+'\n')
save_file.close()
#        #moving the file to another directory
path=os.getcwd()
src = path+"/"+'matches.txt' #source folder
dst = "../res/blast_with_note/match/"#destination folder
shutil.move(src, dst)