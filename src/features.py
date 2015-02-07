from Bio import SeqIO 
from Bio import Entrez
from Bio.Blast import NCBIWWW,NCBIXML
import shutil#moving files
import os.path#cheking files in path
import urllib
from Uniprot_Parser import * #parsing uniprot text file
from Bio.SeqIO import UniprotIO #parsing uniprot xml file
import pandas
import numpy as np


#[0:246000] my zone

#fetching genome zone from Neisseria gonorrhoeae FA 1090 chromosome ncbi
#VERSION: NC_002946.2   GI:59800473


def get_genome_zone(start,stop,filename):
    Entrez.email = "pg27668@alunos.uminho.pt"
    handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", seq_start=start, seq_stop=stop, id="59800473")
    file=open(filename,"w")#creating a GenBank file
    file.write(handle.read())
    file.close()
    handle.close()
    #moving the file to another directory
    path=os.getcwd()
    src = path+"/"+filename #source folder
    dst = "../res"#destination folder
    shutil.move(src, dst)
    record = SeqIO.read("../res/"+filename, "genbank") 
    return record


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


def tabela(lista,locus):
    headers=['locus','gene name','location','GI GeneID','EC','proteinID','Uniprot_ID','product','note']
    dado=[]    
    for i in range(len(locus)):
        dado.append(locus[i])
    data=np.array(lista)
    df=pandas.DataFrame(data, dado, headers)
    df.to_csv("../res/excel/teste", sep='\t')
    
    
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
    

def tabela_pseudogenes(lista,locus):
    headers=['locus','location','gene_ID']
    dado=[]    
    for i in range(len(locus)):
        dado.append(locus[i])
    data=np.array(lista)
    df=pandas.DataFrame(data, dado, headers)
    df.to_csv("../res/excel/teste_pseudogenes", sep='\t')
    

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

#list of genes names
def list_genes_names(record,locus):
    genes=[]
    for i in range(len(locus)):
        if genes_names(record,locus[i])!="Nao tem nome!":
            genes.append(genes_names(record,locus[i]))
    return genes
    
  
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



def tabela_tRNA(lista,locus):
    headers=['locus','gene_ID','product','location']
    dado=[]    
    for i in range(len(locus)):
        dado.append(locus[i])
    data=np.array(lista)
    df=pandas.DataFrame(data, dado, headers)
    df.to_csv("../res/excel/teste_tRNA", sep='\t')
    
              

def uniprot_ID(record):
    locus=locus_tag(record)
    proteins=[]
    IDs=[]
    for i in range(len(locus)):
        proteins.append(protein_ID(record,locus[i]))
    for j in range(len(proteins)):
        data = urllib.request.urlopen("http://www.uniprot.org/uniprot/?query="+proteins[j]+"&sort=score").read()
        if len(data.split())>=3298:
            x=str(data.split()[3298])
            IDs.append(x[6:12])
        else:
            IDs.append(proteins[j])
    return IDs


#['b\'id="Q5FAJ2"\'', 'b\'id="Q5FAJ1"\'', 'b\'id="Q5FAJ0"\'', 'b\'id="Q5FAI9"\'', 'b\'id="Q5FAJ3"\'', 'YP_009115477.1', 'b\'id="Q5FAK0"\'', 'b\'id="Q5FAJ8"\'', 'b\'id="Q5FAJ7"\'', 'b\'id="Q5FAJ6"\'', 'b\'id="Q5FAJ5"\'', 'b\'id="Q5FAJ4"\'', 'b\'id="Q5FAL3"\'', 'b\'id="Q5FAL2"\'', 'b\'id="Q5FAL1"\'', 'b\'id="Q5FAK9"\'', 'b\'id="Q5FAK8"\'', 'b\'id="Q5FAK7"\'', 'b\'id="Q5FAK6"\'', 'b\'id="Q5FAK5"\'', 'b\'id="Q5FAI8"\'', 'b\'id="Q5FAI7"\'', 'b\'id="Q5FAI6"\'', 'b\'id="Q5FAI5"\'', 'b\'id="Q5FAI4"\'', 'b\'id="Q5FAI3"\'', 'b\'id="Q5FAI2"\'', 'b\'id="Q5FAI1"\'', 'b\'id="Q5FAH9"\'', 'b\'id="Q5FAH8"\'', 'b\'id="Q5FAH7"\'', 'b\'id="Q5FAH6"\'', 'b\'id="Q5FAH5"\'', 'b\'id="Q5FAH4"\'', 'b\'id="Q5FAH3"\'', 'b\'id="Q5FAH2"\'', 'b\'id="Q5FAH1"\'', 'b\'id="Q5FAH0"\'', 'b\'id="Q5FAG9"\'', 'b\'id="Q5FAG7"\'', 'b\'id="Q5FAG6"\'', 'b\'id="Q5FAG5"\'', 'b\'id="Q5FAG4"\'', 'b\'id="Q5FAG3"\'', 'b\'id="Q5FAG2"\'', 'b\'id="Q5FAG1"\'', 'b\'id="Q5FAG0"\'', 'b\'id="Q5FAF9"\'', 'b\'id="Q5FAF8"\'', 'b\'id="Q5FAF7"\'', 'YP_008914846.1', 'b\'id="Q5FAF6"\'', 'b\'id="Q5FAF5"\'', 'b\'id="Q5FAF4"\'', 'b\'id="Q5FAF3"\'', 'b\'id="Q5FAF2"\'', 'b\'id="Q5FAF1"\'', 'b\'id="Q5FAE8"\'', 'b\'id="Q5FAE7"\'', 'b\'id="Q5FAE6"\'', 'b\'id="Q5FAE5"\'', 'b\'id="Q5FAE4"\'', 'b\'id="Q5FAE3"\'', 'b\'id="Q5FAE2"\'', 'b\'id="Q5FAE1"\'', 'b\'id="Q5FAE0"\'', 'b\'id="Q5FAD9"\'', 'b\'id="Q5FAD8"\'', 'b\'id="Q5FAD7"\'', 'b\'id="Q5FAD6"\'', 'b\'id="Q5FAD5"\'', 'b\'id="Q5FAD4"\'', 'b\'id="Q5FAD3"\'', 'b\'id="Q5FAD2"\'', 'b\'id="Q5FAD1"\'', 'b\'id="Q5FAD0"\'', 'b\'id="Q5FAC9"\'', 'b\'id="Q5FAC8"\'', 'b\'id="Q5FAC7"\'', 'b\'id="Q5FAC6"\'', 'b\'id="Q5FAC5"\'', 'b\'id="Q5FAC4"\'', 'b\'id="Q5FAC3"\'', 'b\'id="Q5FAC2"\'', 'b\'id="Q5FAC1"\'', 'b\'id="Q5FAC0"\'', 'b\'id="Q5FAB9"\'', 'b\'id="Q5FAB7"\'', 'b\'id="Q5FAB6"\'', 'b\'id="Q5FAB5"\'', 'b\'id="Q5FAB4"\'', 'b\'id="Q5FAB3"\'', 'b\'id="Q5FAB1"\'', 'b\'id="Q5FAB0"\'', 'b\'id="Q5FAA9"\'', 'b\'id="Q5FAA8"\'', 'b\'id="Q5FAA7"\'', 'b\'id="Q5FAA6"\'', 'b\'id="Q5FAA5"\'', 'b\'id="Q5FAA4"\'', 'b\'id="Q5FAA3"\'', 'b\'id="Q5FAA2"\'', 'YP_009115478.1', 'b\'id="Q5FAA1"\'', 'b\'id="Q5FA99"\'', 'b\'id="Q5FA98"\'', 'b\'id="Q5FA97"\'', 'b\'id="Q5FA96"\'', 'b\'id="Q5FA94"\'', 'b\'id="Q5FA93"\'', 'b\'id="Q5FA92"\'', 'b\'id="Q5FA91"\'', 'b\'id="Q5FA90"\'', 'b\'id="Q5FA89"\'', 'b\'id="Q5FA88"\'', 'b\'id="Q5FA87"\'', 'b\'id="Q5FA86"\'', 'b\'id="Q5FA85"\'', 'b\'id="Q5FA84"\'', 'b\'id="Q5FA83"\'', 'YP_008914847.1', 'b\'id="Q5FA80"\'', 'b\'id="Q5FA78"\'', 'YP_207322.2', 'b\'id="Q5FA76"\'', 'b\'id="Q5FA75"\'', 'b\'id="Q5FA73"\'', 'b\'id="Q5FA72"\'', 'b\'id="Q5FA71"\'', 'b\'id="Q5FA70"\'', 'b\'id="Q5FA68"\'', 'b\'id="Q5FA67"\'', 'b\'id="Q5FA66"\'', 'b\'id="Q5FA65"\'', 'b\'id="Q5FA63"\'', 'b\'id="Q5FA62"\'', 'b\'id="Q5FA61"\'', 'b\'id="Q5FA60"\'', 'b\'id="Q5FA59"\'', 'b\'id="Q5FA58"\'', 'b\'id="Q5FA57"\'', 'b\'id="Q5FA56"\'', 'b\'id="Q5FA55"\'', 'b\'id="Q5FA54"\'', 'b\'id="Q5FA53"\'', 'b\'id="Q5FA52"\'', 'b\'id="Q5FA51"\'', 'b\'id="Q5FA50"\'', 'b\'id="Q5FA49"\'', 'b\'id="Q5FA48"\'', 'b\'id="Q5FA47"\'', 'b\'id="Q5FA46"\'', 'b\'id="Q5FA45"\'', 'b\'id="Q5FA44"\'', 'b\'id="Q5FA43"\'', 'b\'id="Q5FA42"\'', 'b\'id="Q5FA41"\'', 'b\'id="Q5FA40"\'', 'b\'id="Q5FA39"\'', 'b\'id="Q5FA38"\'', 'b\'id="Q5FA37"\'', 'b\'id="Q5FA36"\'', 'b\'id="Q5FA35"\'', 'b\'id="Q5FA34"\'', 'b\'id="Q5FA33"\'', 'b\'id="Q5FA32"\'', 'b\'id="Q5FA31"\'', 'b\'id="Q5FA30"\'', 'b\'id="Q5FA29"\'', 'b\'id="Q5FA28"\'', 'b\'id="Q5FA27"\'', 'b\'id="Q5FA26"\'', 'b\'id="Q5FA25"\'', 'b\'id="Q5FA24"\'', 'b\'id="Q5FA23"\'', 'b\'id="Q5FA22"\'', 'b\'id="Q5FA21"\'', 'b\'id="Q5FA20"\'', 'b\'id="Q5FA19"\'', 'b\'id="Q5FA18"\'', 'b\'id="Q5FA17"\'', 'b\'id="Q5FA15"\'', 'b\'id="Q5FA14"\'', 'b\'id="Q5FA13"\'', 'b\'id="Q5FA12"\'', 'b\'id="Q5FA11"\'', 'b\'id="Q5FA10"\'', 'b\'id="Q5FA09"\'', 'b\'id="Q5FA08"\'', 'b\'id="Q5FA07"\'', 'b\'id="Q5FA06"\'', 'b\'id="Q5FA05"\'', 'b\'id="Q5FA04"\'', 'b\'id="Q5FA03"\'', 'b\'id="Q5FA02"\'', 'b\'id="Q5FA01"\'', 'b\'id="Q5FA00"\'', 'b\'id="Q5F9Z9"\'', 'b\'id="Q5F9Z8"\'', 'b\'id="Q5F9Z7"\'', 'b\'id="Q5F9Z5"\'', 'b\'id="Q5F9Z4"\'', 'b\'id="Q5F9Z2"\'']

ID=['Q5FAJ2', 'Q5FAJ1', 'Q5FAJ0', 'Q5FAI9', 'Q5FAJ3', 'YP_009115477.1', 'Q5FAK0', 'Q5FAJ8', 'Q5FAJ7', 'Q5FAJ6', 'Q5FAJ5', 'Q5FAJ4', 'Q5FAL3', 'Q5FAL2', 'Q5FAL1', 'Q5FAK9', 'Q5FAK8', 'Q5FAK7', 'Q5FAK6', 'Q5FAK5', 'Q5FAI8', 'Q5FAI7', 'Q5FAI6', 'Q5FAI5', 'Q5FAI4', 'Q5FAI3', 'Q5FAI2', 'Q5FAI1', 'Q5FAH9', 'Q5FAH8', 'Q5FAH7', 'Q5FAH6', 'Q5FAH5', 'Q5FAH4', 'Q5FAH3', 'Q5FAH2', 'Q5FAH1', 'Q5FAH0', 'Q5FAG9', 'Q5FAG7', 'Q5FAG6', 'Q5FAG5', 'Q5FAG4', 'Q5FAG3', 'Q5FAG2', 'Q5FAG1', 'Q5FAG0', 'Q5FAF9', 'Q5FAF8', 'Q5FAF7', 'YP_008914846.1', 'Q5FAF6', 'Q5FAF5', 'Q5FAF4', 'Q5FAF3', 'Q5FAF2', 'Q5FAF1', 'Q5FAE8', 'Q5FAE7', 'Q5FAE6', 'Q5FAE5', 'Q5FAE4', 'Q5FAE3', 'Q5FAE2', 'Q5FAE1', 'Q5FAE0', 'Q5FAD9', 'Q5FAD8', 'Q5FAD7', 'Q5FAD6', 'Q5FAD5', 'Q5FAD4', 'Q5FAD3', 'Q5FAD2', 'Q5FAD1', 'Q5FAD0', 'Q5FAC9', 'Q5FAC8', 'Q5FAC7', 'Q5FAC6', 'Q5FAC5', 'Q5FAC4', 'Q5FAC3', 'Q5FAC2', 'Q5FAC1', 'Q5FAC0', 'Q5FAB9', 'Q5FAB7', 'Q5FAB6', 'Q5FAB5', 'Q5FAB4', 'Q5FAB3', 'Q5FAB1', 'Q5FAB0', 'Q5FAA9', 'Q5FAA8', 'Q5FAA7', 'Q5FAA6', 'Q5FAA5', 'Q5FAA4', 'Q5FAA3', 'Q5FAA2', 'YP_009115478.1', 'Q5FAA1', 'Q5FA99', 'Q5FA98', 'Q5FA97', 'Q5FA96', 'Q5FA94', 'Q5FA93', 'Q5FA92', 'Q5FA91', 'Q5FA90', 'Q5FA89', 'Q5FA88', 'Q5FA87', 'Q5FA86', 'Q5FA85', 'Q5FA84', 'Q5FA83', 'YP_008914847.1', 'Q5FA80', 'Q5FA78', 'YP_207322.2', 'Q5FA76', 'Q5FA75', 'Q5FA73', 'Q5FA72', 'Q5FA71', 'Q5FA70', 'Q5FA68', 'Q5FA67', 'Q5FA66', 'Q5FA65', 'Q5FA63', 'Q5FA62', 'Q5FA61', 'Q5FA60', 'Q5FA59', 'Q5FA58', 'Q5FA57', 'Q5FA56', 'Q5FA55', 'Q5FA54', 'Q5FA53', 'Q5FA52', 'Q5FA51', 'Q5FA50', 'Q5FA49', 'Q5FA48', 'Q5FA47', 'Q5FA46', 'Q5FA45', 'Q5FA44', 'Q5FA43', 'Q5FA42', 'Q5FA41', 'Q5FA40', 'Q5FA39', 'Q5FA38', 'Q5FA37', 'Q5FA36', 'Q5FA35', 'Q5FA34', 'Q5FA33', 'Q5FA32', 'Q5FA31', 'Q5FA30', 'Q5FA29', 'Q5FA28', 'Q5FA27', 'Q5FA26', 'Q5FA25', 'Q5FA24', 'Q5FA23', 'Q5FA22', 'Q5FA21', 'Q5FA20', 'Q5FA19', 'Q5FA18', 'Q5FA17', 'Q5FA15', 'Q5FA14', 'Q5FA13', 'Q5FA12', 'Q5FA11', 'Q5FA10', 'Q5FA09', 'Q5FA08', 'Q5FA07', 'Q5FA06', 'Q5FA05', 'Q5FA04', 'Q5FA03', 'Q5FA02', 'Q5FA01', 'Q5FA00', 'Q5F9Z9', 'Q5F9Z8', 'Q5F9Z7', 'Q5F9Z5', 'Q5F9Z4', 'Q5F9Z2']


def info_uniprot():
    identifier=[]
    handle = open("../res/"+"uniprot.txt")
    records = parse(handle) # Uses the function 'parse' from the module. 
    for record in records:
        for i in range(len(ID)):
            identifier.append([])
            if record["AC"]==ID[i]+";":
                #identifier.append([])
                identifier[i].append(ID[i])
                identifier[i].append(record["ID"])
                identifier[i].append(record["DE"])
                identifier[i].append(record["CC"])
    handle.close()
    return identifier


def more_info_uniprot():
    handle = open("../res/"+"uniprot_XML.xml")
    records=UniprotIO.UniprotIterator(handle,return_raw_comments=True)
    refs=[]
    for record in records:
        for i in range(len(ID)):
            if record.id==ID[i]:
                refs.append([ID[i]]+record.dbxrefs)#GO´s
    handle.close()
    return refs
    

#creating data frame for posterior excel file
def tabela_uniprot():
    dado=[]
    lista=[]
    ident=info_uniprot()
    for i in range(len(ident)):
        if ident[i]!=[]:
            dado.append(ident[i][0])       
            lista.append(ident[i])
    data=np.array(lista)
    df=pandas.DataFrame(data, dado)
    df.to_csv("../res/excel/teste_uniprot", sep='\t')
    

#creating data frame for posterior excel file   
def tabela_uniprot2():
    dado=[]
    lista=[]
    refs=sorting(tab())
    for i in range(len(refs)):
        dado.append(refs[i][0])
        lista.append(refs[i])
    data=np.array(lista)
    df=pandas.DataFrame(data, dado)
    df.to_csv("../res/excel/teste_uniprot2", sep='\t')


#appending extra columns for posterior data frame
def tab():
    refs=more_info_uniprot()
    lista=[]
    l=[]
    for i in range(len(refs)):
        l.append(len(refs[i]))
        lista.append([])
    m=max(l)
    for j in range(len(refs)):
        if len(refs[j])<m:
            for k in range(len(refs[j])):
                lista[j].append(refs[j][k])
            for h in range(len(refs[j]),m):
                lista[j].append("nao tem")
        elif len(refs[j])==m:
            for g in range(len(refs[j])):
                lista[j].append(refs[j][g])        
    return lista
    
    
#sorting by protein ID
def sorting(lista):
    mat=[]
    for i in range(len(ID)):
        for j in range(len(lista)):
            if ID[i]==lista[j][0]:
                mat.append(lista[j])
    return mat


def GO():
    l=[]
    lista=sorting(tab())
    for i in range(len(lista)):
        for j in range(len(lista[i])):
            if "GO:GO:" in lista[i][j]:
                l.append(lista[i][0])
                l.append(lista[i][j])
    return l


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


#745998704  
#Running Blast and saving info into a file
#GI_number - gene identification number
def blast(GI_numb,filename):
    result_handle = NCBIWWW.qblast("blastp", "swissprot", GI_numb)
    save_file = open(filename, "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()
    #moving the file to another directory
    path=os.getcwd()
    src = path+"/"+filename #source folder
    dst = "../res"#destination folder
    shutil.move(src, dst)
    
    
#Parsing Blast files
def parse_blast(filename):
    #E_VALUE_THRESH = 0.05
    result_handle = open("../res/"+filename)
    blast_record = NCBIXML.read(result_handle)
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            #if hsp.expect < E_VALUE_THRESH:
            print ("****Alignment****")
            print ('sequence:', alignment.title)
            print ('length:', alignment.length)
            print ('e value:', hsp.expect)
    result_handle.close()


#Get gi from protein without note    
def giwithout_note(record):
    ID=[]
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
            if "note" not in my_cds.qualifiers:
                x=my_cds.qualifiers["db_xref"]
                gi=x[0]
                ID.append(gi[3:])
    return ID

#blast gi without note    
def blastnote(filename):
    gi=giwithout_note(record)
    for i in range(44,len(gi)):
        GI_numb=str(gi[i])
        result_handle = NCBIWWW.qblast("blastp","swissprot", GI_numb)
        save_file = open(GI_numb+'.xml', "w")
        save_file.write(result_handle.read())
        save_file.close()
        result_handle.close()
        #moving the file to another directory
        path=os.getcwd()
        src = path+"/"+GI_numb+'.xml' #source folder
        dst = "../res/blast_without_note"#destination folder
        shutil.move(src, dst)
        
        
def blastanaliser():
    blast=[]    
    for file in os.listdir("../res/blast_without_note"):
        if file.endswith(".xml"):
            blast.append(file)
    E_VALUE_THRESH = 0.05
    lista=[]
    for i in range(len(blast)):
        lista.append([])
        lista[i].append(blast[i])
        result_handle = open("../res/blast_without_note/"+blast[i])
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
    dst = "../res/blast_without_note/nomatch/"#destination folder
    shutil.move(src, dst)
                     
    save_file = open('matches.txt', "w")
    for i in range(len(lista)):                
        if len(lista[i])>2:
            save_file.write(str(lista[i])+'\n')
    save_file.close()
    #        #moving the file to another directory
    path=os.getcwd()
    src = path+"/"+'matches.txt' #source folder
    dst = "../res/blast_without_note/match/"#destination folder
    shutil.move(src, dst)

def gimatch():
    handle = open("../res/blast_without_note/match/matches.txt").readlines()    
    gimatch=[]
    gi=[]
    hit=[]
    s2='sp|'
    xml='['
    for i in handle:
        n = i[i.index(xml) + len(xml):] 
        g = n.split('.xml', 1)[0]  
        gi.append(g[1:])
    for J in handle:
        s3 = J[J.index(s2) + len(s2):] 
        sP = s3.split('|', 1)[0] 
        sp=sP[:6]    
        hit.append(sp)
    for k  in range(len(gi)):
        gimatch.append(gi[k]+' '+hit[k])    
            
            
    
    for n in range(len(gi)):
        site = urllib.request.urlopen("http://www.uniprot.org/uniprot/"+hit[n]+".txt")
        data = site.readlines()
        nome=gi[n]
        file = open("../res/blast_without_note/match/function/"+nome+'.txt',"wb") #open file in binary mode
        file.writelines(data)
        file.close()
        
        
#function to see possible function of proteins that didnt have note:
def getfunction():
    blast=[]    
    for file in os.listdir("../res/blast_without_note/match/function/"):
        if file.endswith(".txt"):
               blast.append(file)
    
    lista=[]
    for j in range(len(blast)):
        lista.append([])
        nome=blast[j]
        gi = nome.replace(".txt","")
        first = '-!- FUNCTION:'
        last = 'CC'
        file = open("../res/blast_without_note/match/function/"+nome).read()
        data = file.replace("\n", " ") 
        try:
            start = data.rindex( first ) + len( first )
            end = data.rindex(last, start)
            novo= data[start:end] 
            lista[j].append('Gi: '+gi+'  '+'Possivel função:  '+novo)
        except:
            pass   
    return lista
    
#return all hits from blast
def allhits():
    lista=[]
    handle = open("../res/blast_without_note/match/matches.txt").readlines()    
    first = 'sp|'
    last = '|'
    gi=[]
    xml='['
    
    
        
    for q in range(len(handle)):
        lista.append([])
        x=handle[q].split()
        for i in range (len(x)):
            m=x[i]
            if 'sp|' in m:    
                try:
                    start = m.rindex( first ) + len( first )
                    end = m.rindex( last, start )
                    novo= m[start:end]
                    lista[q].append(novo)
                except:
                    pass
          
    
    for f in handle:
        n = f[f.index(xml) + len(xml):] 
        g = n.split('.xml', 1)[0]  
        gi.append(g[1:])
    
    for p in range(len(lista)):
        lista[p].append(gi[p])
    return lista

def uniprotallhits():
    handle = open("../res/blast_without_note/match/allhits/allhits.txt").readlines()
    
    
    
    for n in range(len (handle)):       
            x=handle[n].split()
            for k in range(len(x)-1,len(x)):
                e=x[k]
                limpo=e.replace("']","")
                gi=(limpo[1:])
                if not os.path.exists("../res/blast_without_note/match/function/teste/"+gi):
                    os.makedirs("../res/blast_without_note/match/function/teste/"+gi)
            for j in range(len(x)-1):
                m=x[j]
                q=m.replace("[" ,"")
                protein=q[1:7]
                site = urllib.request.urlopen("http://www.uniprot.org/uniprot/"+protein+".txt")
                data = site.readlines()
                file = open("../res/blast_without_note/match/function/teste/"+protein+'.txt',"wb") #open file in binary mode
                file.writelines(data)
                file.close()
                try:
                    src = "../res/blast_without_note/match/function/teste/"+protein+'.txt' #source folder
                    dst = "../res/blast_without_note/match/function/teste/"+gi #destination folder
                    shutil.move(src, dst)
                except:
                        pass
def allfunctions():
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
                        
                    file = open("../res/blast_without_note/match/allhits/funcao_all_hits/"+ginote+".txt",'w')       
                    for i in range(len(lista)):
                        file.write("%s\n" % lista[i])
                    file.close()
                
    
def menu(record):
    ans=True
    while ans:
        print("""
    1.Accessing info from one gene
    2.without note (returns proteinID)
    3.translation (needs proteinID)
    4.tRNA
    5.pseudogenes
    6.hypothetical proteins
    7.Article with a gene reference
    8.Running protein Blast (needs GI number)
    9.Parsing blast
    10.Uniprot_ID 
    11.Uniprot Identifier, Definition, Subcellular location
    12.Get gi from protein without note
    13.tabela
    14.Get note from protein without note
    15.blastanaliser
    16.pega no primeiro hit do blast
    17.GO numbers
    18.function to see possible function of proteins that didnt have note:
    19.return all hits from blast
    20.Go to uniprot and download information for all hits
    21.Get all information for all hits
    60.Exit
    """)
        ans=input("Choose an option? ")
        if ans=="1":
            locus=locus_tag(record)            
            fa=info(record,locus)
            nr=int(input("Qual gene?"))
            print(aceder(fa,nr))
#            print(locus_tag(record))
#            print(genes_names(record))
        elif ans=="2":
            print(without_note(record))
        elif ans=="3":
            prot_ID=str(input("Protein ID: "))
            print(translation(record,prot_ID))
        elif ans=="4":
            loc=tRNA(record)
            print(info_tRNA(record,loc))
        elif ans=="5":
            locus=pseudogenes(record)
            print(info_pseudogenes(record,locus))
        elif ans=="6":
            print(hypoth_proteins(record))
        elif ans=="7":
            gene=str(input("Gene name: "))
            print(DB_pubmed(gene))
        elif ans=="8":
            file=str(input("Qual o nome a colocar no ficheiro? "))+".xml"
            while os.path.isfile("../res/"+file):
                 file=str(input("Qual o nome a colocar no ficheiro? "))+".xml"
            else:False
            GI=str(input("Qual o GI number da sequencia? "))
            blast(GI,file)
        elif ans=="9":
            file=str(input("Qual o nome do ficheiro? "))+".xml"
            parse_blast(file)
        elif ans=="10":
            #print(uniprot_ID(record))
            print(ID)
        elif ans=="11":
            print(info_uniprot())
            print(more_info_uniprot())
        elif ans=="12":
            print(giwithout_note(record))
        elif ans=="13":
            locus=locus_tag(record)            
            fa=info(record,locus)
            print(tabela(fa,locus))
            locus_tRNA=tRNA(record)
            lista_tRNA=info_tRNA(record,locus_tRNA)
            print(tabela_tRNA(lista_tRNA,locus_tRNA))
            locus_pseudo=pseudogenes(record)
            lista_pseudo=info_pseudogenes(record,locus_pseudo)
            print(tabela_pseudogenes(lista_pseudo,locus_pseudo))
            print(tabela_uniprot())
            print(tabela_uniprot2())
        elif ans=="14":
            blastnote(filename)
        elif ans=="15":
            blastanaliser()
        elif ans=="16":
            gimatch()
        elif ans=="17":
            print(GO())
            
        elif ans=="18":
           file = open("../res/blast_without_note/match/allhits/funcao.txt",'w')
           lista=getfunction()        
           for i in range(len(lista)):
               file.write("%s\n" % lista[i])
#           print(getfunction())
        elif ans=="19":
           file = open("../res/blast_without_note/match/allhits/allhits.txt",'w')
           lista=allhits()         
           for i in range(len(lista)):
               file.write("%s\n" % lista[i])
           file.close()
        elif ans=="20":
            uniprotallhits()
            
        elif ans=="21":    
            allfunctions()
        elif ans=="60":
            ans = False
        else:
            print("\nInvalid")

#gimatch()

def create_file():
    start=str(input(print("Insira o inicio da sua zona do genoma:")))
    stop=str(input(print("Insira o fim da sua zona do genoma:")))
    filename=str(input(print("Insira nome do ficheiro: "))+".gb")#filename plus extension being genbank
    #checking if name already exists    
    while os.path.isfile("../res/"+filename):
        filename=str(input(print("Insira outro nome: "))+".gb")
    else:False
    #if not exists get the genome zone and create new file
    record=get_genome_zone(start,stop,filename)
    return record
    


#main
if __name__ == "__main__":
    res=str(input(print("Ja tem ficheiro S/N? "))).upper()
    if res=="N":
        record=create_file()
    elif res=="S":
        filename=str(input("Nome do ficheiro? "))+".gb"
        record = SeqIO.read("../res/"+filename, "genbank") 
    menu(record)
    
