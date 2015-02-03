# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 23:45:12 2015

@author: Danielbraga
"""
##FUNCTION: 
##def function():
#handle = open("uniprot.txt").readlines()
#function='CC   '
#uniprotid='AC' 
#listfunction=[]
#uniprot=[]
#lista=[]
#
##for i in handle:
##    if function in i:
##        y=i.split(" ")             
##        listfunction.append(i)
##
##for j in handle:
##    if uniprotid in j:
##        uniprot.append(j)
##x=len(listname)        
##for k in range(x):
##    lista.append(((uniprot[k]+" "+listname[k])))
##return lista
###print(review())
#
#
from Bio import SeqIO 
from Bio import Entrez
from Bio.Blast import NCBIWWW,NCBIXML
import shutil#moving files
import os.path#cheking files in path
import urllib
from Uniprot_Parser import * #parsing uniprot text file
from Bio.SeqIO import UniprotIO #parsing uniprot xml file

lista=['YP_207189.1', 'YP_009115477.1', 'YP_207195.1', 'YP_207197.1', 'YP_207199.1', 'YP_207200.1', 'YP_207201.1', 'YP_207202.1', 'YP_207203.1', 'YP_207204.1', 'YP_207206.1', 'YP_207207.1', 'YP_207209.1', 'YP_207211.1', 'YP_207212.1', 'YP_207213.1', 'YP_207214.1', 'YP_207216.1', 'YP_207223.1', 'YP_207225.1', 'YP_207229.1', 'YP_207233.1', 'YP_207234.1', 'YP_207235.1', 'YP_207236.1', 'YP_207237.1', 'YP_207240.1', 'YP_207241.1', 'YP_207242.1', 'YP_008914846.1', 'YP_207245.1', 'YP_207248.1', 'YP_207251.1', 'YP_207252.1', 'YP_207253.1', 'YP_207254.1', 'YP_207257.1', 'YP_207260.1', 'YP_207261.1', 'YP_207262.1', 'YP_207263.1', 'YP_207264.1', 'YP_207267.1', 'YP_207268.1', 'YP_207269.1', 'YP_207271.1', 'YP_207274.1', 'YP_207275.1', 'YP_207276.1', 'YP_207277.1', 'YP_207282.1', 'YP_207283.1', 'YP_207284.1', 'YP_207286.1', 'YP_207291.1', 'YP_207292.1', 'YP_207293.1', 'YP_207294.1', 'YP_207295.1', 'YP_207296.1', 'YP_207297.1', 'YP_009115478.1', 'YP_207301.1', 'YP_207302.1', 'YP_207303.1', 'YP_207306.1', 'YP_207307.1', 'YP_207308.1', 'YP_207309.1', 'YP_207310.1', 'YP_207311.1', 'YP_207312.1', 'YP_207313.1', 'YP_207314.1', 'YP_207315.1', 'YP_008914847.1', 'YP_207319.1', 'YP_207321.1', 'YP_207322.2', 'YP_207326.1', 'YP_207327.1', 'YP_207328.1', 'YP_207329.1', 'YP_207331.1', 'YP_207332.1', 'YP_207333.1', 'YP_207334.1', 'YP_207337.1', 'YP_207338.1', 'YP_207343.1', 'YP_207344.1', 'YP_207345.1', 'YP_207346.1', 'YP_207348.1', 'YP_207349.1', 'YP_207353.1', 'YP_207354.1', 'YP_207355.1', 'YP_207359.1', 'YP_207360.1', 'YP_207361.1', 'YP_207363.1', 'YP_207366.1', 'YP_207367.1', 'YP_207368.1', 'YP_207369.1', 'YP_207370.1', 'YP_207371.1', 'YP_207373.1', 'YP_207379.1', 'YP_207380.1', 'YP_207381.1', 'YP_207385.1', 'YP_207387.1', 'YP_207391.1', 'YP_207392.1', 'YP_207393.1', 'YP_207395.1', 'YP_207397.1', 'YP_207398.1', 'YP_207400.1', 'YP_207401.1', 'YP_207402.1']
record = SeqIO.read("../res/"+filename, "genbank") 
def gene_ID_GI(record,lista):
    for i in range(len(record.features)):
        my_gene = record.features[i]
        if my_gene.type == "CDS":             
        #there are usually several db_xref entries in CDS, but one in gene qualifier
            if my_gene.qualifiers["locus_tag"][0]==str(locus_tag):    
                if "db_xref" in my_gene.qualifiers:
                    x=my_gene.qualifiers["db_xref"]
                    return x[0],x[1]
        
