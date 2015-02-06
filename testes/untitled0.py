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
for file in os.listdir("../res/blast_without_note/match/function"):
    if file.endswith(".txt"):
        blast.append(file)