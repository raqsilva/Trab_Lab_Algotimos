# -*- coding: utf-8 -*-
from Bio import SeqIO #reading gb file 


#Getting pseudogenes locus_tag
def pseudogenes(record):
    locus=[]
    for i in range(len(record.features)):
        my_gene = record.features[i]
        if my_gene.type == "gene":
             if "pseudogene" in my_gene.qualifiers:
                locus.append(my_gene.qualifiers["locus_tag"][0])
    return locus


def location_pseudo(record,locus_tag):
    featcds = [ ]
    for i in range(len(record.features)):
        if record.features[i].type == "gene":
            if record.features[i].qualifiers["locus_tag"][0]==str(locus_tag): 
                featcds.append(i)
    for k in featcds: 
        return record.features[k].location


def pseudogeneID(record,locus_tag):
    for i in range(len(record.features)):
        my_gene = record.features[i]
        if my_gene.type == "gene":
            if my_gene.qualifiers["locus_tag"][0]==str(locus_tag): 
                if "pseudogene" in my_gene.qualifiers:
                    return my_gene.qualifiers["db_xref"][0]


def info_pseudogenes(record,locus):
    lista=[]
    for i in range(len(locus)):
        lista.append([])
        lista[i].append(locus[i])
        lista[i].append(location_pseudo(record,locus[i]))
        lista[i].append(pseudogeneID(record,locus[i]))
    return lista   
    
