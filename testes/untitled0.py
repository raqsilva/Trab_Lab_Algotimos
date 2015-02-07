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

handle = open("../res/blast_without_note/match/allhits/allhits.txt").readlines()



for n in range(len (handle)):       
        x=handle[n].split()
        for k in range(len(x)-1,len(x)):
            e=x[k]
            limpo=e.replace("']","")
            ginote=(limpo[1:])
            blast=[]    
            for file in os.listdir("../res/blast_without_note/match/function/teste/"+ginote):
                if file.endswith(".txt"):
                       blast.append(file)
                lista=[]
                for j in range(len(blast)):
                    nome=blast[j]
                    gi = nome.replace(".txt","")
                    print(ginote)
                    first = '-!- FUNCTION:'
                    last = 'CC'
                    file = open("../res/blast_without_note/match/function/teste/"+ginote+'/'+gi+'.txt').read()
                    data = file.replace("\n", " ") 
                    try:
                        start = data.rindex( first ) + len( first )
                        end = data.rindex(last, start)
                        novo= data[start:end] 
                        lista.append('Gi: '+gi+'  '+'Possivel função:  '+novo)
                    except:
                        pass   
                    
                file = open("../res/blast_without_note/match/function/teste/"+ginote+".txt",'w')       
                for i in range(len(lista)):
                    file.write("%s\n" % lista[i])
                file.close()
                