# -*- coding: utf-8 -*-
from Bio import SeqIO #reading gb file 

#locus from tRNA             
def tRNA(record):
    #trna=[]
    locus=[]
    for i in range(len(record.features)):
        my_trna = record.features[i]
        if my_trna.type == "tRNA":
            #trna.append(my_trna.qualifiers["product"][0])
            locus.append(my_trna.qualifiers["locus_tag"][0])
    return locus


#gets geneID form tRNA
def gene_ID_tRNA(record,locus_tag):
    for i in range(len(record.features)):
        my_gene = record.features[i]
        if my_gene.type == "tRNA":             
            if my_gene.qualifiers["locus_tag"][0]==str(locus_tag):    
                x=my_gene.qualifiers["db_xref"]
                return x[0]
                    

#gets product from tRNA
def product_tRNA(record,locus_tag):
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "tRNA":
            if my_cds.qualifiers["locus_tag"][0]==str(locus_tag):  
                return my_cds.qualifiers["product"][0]
                    

#Gets location strand from tRNA
def location_tRNA(record,locus_tag):
    featcds = [ ]
    for i in range(len(record.features)):
        if record.features[i].type == "tRNA":
            if record.features[i].qualifiers["locus_tag"][0]==str(locus_tag): 
                featcds.append(i)
    for k in featcds: 
        return record.features[k].location


#Gets all the info from tRNA                
def info_tRNA(record,locus):
    lista=[]
    for i in range(len(locus)):
        lista.append([])
        lista[i].append(locus[i])
        lista[i].append(gene_ID_tRNA(record,locus[i]))
        lista[i].append(product_tRNA(record,locus[i]))
        lista[i].append(location_tRNA(record,locus[i]))
    return lista   
