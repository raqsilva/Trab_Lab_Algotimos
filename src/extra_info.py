# -*- coding: utf-8 -*-
from Bio import SeqIO #reading gb file 
from gb_features import genes_names
from Bio import Entrez #fetching genbank file


#proteins without notes
def without_note(record):
    ID=[]
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
            if "note" not in my_cds.qualifiers:
                ID.append(my_cds.qualifiers["protein_id"][0])
    return ID


#get hypothetical proteins (proteinID) from my zone
def hypoth_proteins(record):
    hypo=[]
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
            if my_cds.qualifiers["product"]==["hypothetical protein"]:
                hypo.append(my_cds.qualifiers["protein_id"][0])
    return hypo
   

def GI_number(record,locus_tag):
    GI=[]
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
            for j in range(len(locus_tag)):
                if my_cds.qualifiers["locus_tag"][0]==str(locus_tag[j]):
                    if "db_xref" in my_cds.qualifiers:
                        x=my_cds.qualifiers["db_xref"]
                        GI.append(x[0])
    return GI


#list of genes names
def list_genes_names(record,locus):
    genes=[]
    for i in range(len(locus)):
        if genes_names(record,locus[i])!="Nao tem nome!":
            genes.append(genes_names(record,locus[i]))
    return genes
  
  
#Translated sequence, protein sequence
def translation(record,protein_ID): 
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
             if my_cds.qualifiers["protein_id"][0]==str(protein_ID):
                    if "translation" in my_cds.qualifiers:
                        print (my_cds.qualifiers["translation"][0])
                    else:
                        print( "Nao ha traducao")


#Searching articles from PubMed DB referring to my organism and a gene
def DB_pubmed(gene):
    """
    handle = Entrez.egquery(term="Neisseria gonorrhoeae")
    record = Entrez.read(handle)
    for row in record["eGQueryResult"]:
        if row["DbName"]=="pubmed":
            print(row["Count"])
    """
    handle = Entrez.esearch(db="pubmed", term="Neisseria gonorrhoeae[Orgn] AND "+gene+"[Gene]", retmax=11117)
    record = Entrez.read(handle)
    idlist = record["IdList"]
    return idlist


          