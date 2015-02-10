# -*- coding: utf-8 -*-
from Bio import SeqIO #reading gb file 


#[0:246000] my zone

#fetching genome zone from Neisseria gonorrhoeae FA 1090 chromosome ncbi
#VERSION: NC_002946.2   GI:59800473

#Get all the info from one gene
def info(record,locus):
    lista=[]
    for i in range(len(locus)):
        lista.append([])
        lista[i].append(locus[i])
        lista[i].append(genes_names(record,locus[i]))
        lista[i].append(location(record,locus[i]))
        lista[i].append(gene_ID_GI(record,locus[i]))
        lista[i].append(EC_number(record,locus[i]))
        lista[i].append(protein_ID(record,locus[i]))
        lista[i].append(ID[i])
        lista[i].append(product(record,locus[i]))
        lista[i].append(note(record,locus[i]))
    return lista
    
    
def aceder(lista,nr):
    return lista[nr-1]


#return protein ID
def protein_ID(record,locus_tag):
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
            if my_cds.qualifiers["locus_tag"][0]==str(locus_tag):
                if "protein_id" in my_cds.qualifiers:
                    return my_cds.qualifiers["protein_id"][0]
                else:
                    return "Nao contem protein_id"


#Gives the sequence location in genome
def location(record,locus_tag):
    featcds = [ ]
    for i in range(len(record.features)):
        if record.features[i].type == "CDS":
            if record.features[i].qualifiers["locus_tag"][0]==str(locus_tag): 
                featcds.append(i)
    for k in featcds: 
        return record.features[k].location
                
        
#Notes from the sequence
def note(record,locus_tag):
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
            if my_cds.qualifiers["locus_tag"][0]==str(locus_tag):  
                if "note" in my_cds.qualifiers:
                    return my_cds.qualifiers["note"][0]
                else:
                    return "Nao contem nota!"


#proteins without notes
def without_note(record):
    ID=[]
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
            if "note" not in my_cds.qualifiers:
                ID.append(my_cds.qualifiers["protein_id"][0])
    return ID


#locus tag, genes
def locus_tag(record):
    locus=[]
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS": 
            locus.append(my_cds.qualifiers["locus_tag"][0])
    return locus
    
    
#Products from the sequence, name of the proteins
def product(record,locus_tag):
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
            if my_cds.qualifiers["locus_tag"][0]==str(locus_tag):  
                if "product" in my_cds.qualifiers:
                    return my_cds.qualifiers["product"][0]
                else:
                    return "Nao contem produtos!"


#get hypothetical proteins (proteinID) from my zone
def hypoth_proteins(record):
    hypo=[]
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
            if my_cds.qualifiers["product"]==["hypothetical protein"]:
                hypo.append(my_cds.qualifiers["protein_id"][0])
    return hypo
    
       
#YP_009115478.1  
#Gene ID and GI number
def gene_ID_GI(record,locus_tag):
    for i in range(len(record.features)):
        my_gene = record.features[i]
        if my_gene.type == "CDS":             
        #there are usually several db_xref entries in CDS, but one in gene qualifier
            if my_gene.qualifiers["locus_tag"][0]==str(locus_tag):    
                if "db_xref" in my_gene.qualifiers:
                    x=my_gene.qualifiers["db_xref"]
                    return x[0],x[1]


#Getting genes names
def genes_names(record,locus_tag):
    for i in range(len(record.features)):
        my_gene = record.features[i]
        if my_gene.type == "CDS":
            if my_gene.qualifiers["locus_tag"][0]==str(locus_tag):  
                if "gene" in my_gene.qualifiers:
                    return my_gene.qualifiers["gene"][0]
                else:
                    return "Nao tem nome!"


#protein EC number, identification
def EC_number(record,locus_tag):
    for i in range(len(record.features)):
        my_gene = record.features[i]
        if my_gene.type == "CDS":
            if my_gene.qualifiers["locus_tag"][0]==str(locus_tag):  
                if "EC_number" in my_gene.qualifiers:
                    return my_gene.qualifiers["EC_number"][0]
                else:
                    return "Nao contem EC_number!"


ID=['Q5FAJ2', 'Q5FAJ1', 'Q5FAJ0', 'Q5FAI9', 'Q5FAJ3', 'YP_009115477.1', 'Q5FAK0', 'Q5FAJ8', 'Q5FAJ7', 'Q5FAJ6', 'Q5FAJ5', 'Q5FAJ4', 'Q5FAL3', 'Q5FAL2', 'Q5FAL1', 'Q5FAK9', 'Q5FAK8', 'Q5FAK7', 'Q5FAK6', 'Q5FAK5', 'Q5FAI8', 'Q5FAI7', 'Q5FAI6', 'Q5FAI5', 'Q5FAI4', 'Q5FAI3', 'Q5FAI2', 'Q5FAI1', 'Q5FAH9', 'Q5FAH8', 'Q5FAH7', 'Q5FAH6', 'Q5FAH5', 'Q5FAH4', 'Q5FAH3', 'Q5FAH2', 'Q5FAH1', 'Q5FAH0', 'Q5FAG9', 'Q5FAG7', 'Q5FAG6', 'Q5FAG5', 'Q5FAG4', 'Q5FAG3', 'Q5FAG2', 'Q5FAG1', 'Q5FAG0', 'Q5FAF9', 'Q5FAF8', 'Q5FAF7', 'YP_008914846.1', 'Q5FAF6', 'Q5FAF5', 'Q5FAF4', 'Q5FAF3', 'Q5FAF2', 'Q5FAF1', 'Q5FAE8', 'Q5FAE7', 'Q5FAE6', 'Q5FAE5', 'Q5FAE4', 'Q5FAE3', 'Q5FAE2', 'Q5FAE1', 'Q5FAE0', 'Q5FAD9', 'Q5FAD8', 'Q5FAD7', 'Q5FAD6', 'Q5FAD5', 'Q5FAD4', 'Q5FAD3', 'Q5FAD2', 'Q5FAD1', 'Q5FAD0', 'Q5FAC9', 'Q5FAC8', 'Q5FAC7', 'Q5FAC6', 'Q5FAC5', 'Q5FAC4', 'Q5FAC3', 'Q5FAC2', 'Q5FAC1', 'Q5FAC0', 'Q5FAB9', 'Q5FAB7', 'Q5FAB6', 'Q5FAB5', 'Q5FAB4', 'Q5FAB3', 'Q5FAB1', 'Q5FAB0', 'Q5FAA9', 'Q5FAA8', 'Q5FAA7', 'Q5FAA6', 'Q5FAA5', 'Q5FAA4', 'Q5FAA3', 'Q5FAA2', 'YP_009115478.1', 'Q5FAA1', 'Q5FA99', 'Q5FA98', 'Q5FA97', 'Q5FA96', 'Q5FA94', 'Q5FA93', 'Q5FA92', 'Q5FA91', 'Q5FA90', 'Q5FA89', 'Q5FA88', 'Q5FA87', 'Q5FA86', 'Q5FA85', 'Q5FA84', 'Q5FA83', 'YP_008914847.1', 'Q5FA80', 'Q5FA78', 'YP_207322.2', 'Q5FA76', 'Q5FA75', 'Q5FA73', 'Q5FA72', 'Q5FA71', 'Q5FA70', 'Q5FA68', 'Q5FA67', 'Q5FA66', 'Q5FA65', 'Q5FA63', 'Q5FA62', 'Q5FA61', 'Q5FA60', 'Q5FA59', 'Q5FA58', 'Q5FA57', 'Q5FA56', 'Q5FA55', 'Q5FA54', 'Q5FA53', 'Q5FA52', 'Q5FA51', 'Q5FA50', 'Q5FA49', 'Q5FA48', 'Q5FA47', 'Q5FA46', 'Q5FA45', 'Q5FA44', 'Q5FA43', 'Q5FA42', 'Q5FA41', 'Q5FA40', 'Q5FA39', 'Q5FA38', 'Q5FA37', 'Q5FA36', 'Q5FA35', 'Q5FA34', 'Q5FA33', 'Q5FA32', 'Q5FA31', 'Q5FA30', 'Q5FA29', 'Q5FA28', 'Q5FA27', 'Q5FA26', 'Q5FA25', 'Q5FA24', 'Q5FA23', 'Q5FA22', 'Q5FA21', 'Q5FA20', 'Q5FA19', 'Q5FA18', 'Q5FA17', 'Q5FA15', 'Q5FA14', 'Q5FA13', 'Q5FA12', 'Q5FA11', 'Q5FA10', 'Q5FA09', 'Q5FA08', 'Q5FA07', 'Q5FA06', 'Q5FA05', 'Q5FA04', 'Q5FA03', 'Q5FA02', 'Q5FA01', 'Q5FA00', 'Q5F9Z9', 'Q5F9Z8', 'Q5F9Z7', 'Q5F9Z5', 'Q5F9Z4', 'Q5F9Z2']

