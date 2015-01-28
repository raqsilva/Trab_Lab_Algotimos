from Bio import SeqIO 
from Bio import Entrez
from Bio.Blast import NCBIWWW,NCBIXML
import shutil#moving files
import os.path#cheking files in path

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
        lista[i].append(note(record,locus[i]))
        lista[i].append(product(record,locus[i]))
        lista[i].append(gene_ID_GI(record,locus[i]))
        lista[i].append(EC_number(record,locus[i]))
    return lista
    
    
def aceder(lista,nr):
    return lista[nr]


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
        

#Getting pseudogenes
def pseudogenes(record):
    genes=[]
    for i in range(len(record.features)):
        my_gene = record.features[i]
        if my_gene.type == "gene": 
            if "pseudogene" in my_gene.qualifiers:
                genes.append(my_gene.qualifiers["db_xref"][0])
    return genes


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
                    

#information from tRNA             
def tRNA(record):
    trna=[]
    locus=[]
    for i in range(len(record.features)):
        my_trna = record.features[i]
        if my_trna.type == "tRNA":
            trna.append(my_trna.qualifiers["product"][0])
            locus.append(my_trna.qualifiers["locus_tag"][0])
    return trna,locus
  

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



def menu(record):
    ans=True
    while ans:
        print("""
    1.locus_tag plus genes names
    2.product
    3.note (returns note) and without note (returns proteinID)
    4.GI number, geneID (needs proteinID)
    5.translation (needs proteinID)
    6.EC_number
    7.location
    8.tRNA
    9.pseudogenes
    10.hypothetical proteins
    11.Article with a gene reference
    12.Running protein Blast (needs GI number)
    13.Parsing blast
    14.Exit
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
            print(product(record))
        elif ans=="3":
            print(without_note(record))
            print(note(record))
        elif ans=="4":
            prot_ID=str(input("Protein ID: "))
            gene_ID_GI(record,prot_ID)
        elif ans=="5":
            prot_ID=str(input("Protein ID: "))
            print(translation(record,prot_ID))
        elif ans=="6":
            print(EC_number(record))
        elif ans=="7":
            location(record)
        elif ans=="8":
            print(tRNA(record))
        elif ans=="9":
            print(pseudogenes(record))
        elif ans=="10":
            print(hypoth_proteins(record))
        elif ans=="11":
            gene=str(input("Gene name: "))
            print(DB_pubmed(gene))
        elif ans=="12":
            file=str(input("Qual o nome a colocar no ficheiro? "))+".xml"
            while os.path.isfile("../res/"+file):
                 file=str(input("Qual o nome a colocar no ficheiro? "))+".xml"
            else:False
            GI=str(input("Qual o GI number da sequencia? "))
            blast(GI,file)
        elif ans=="13":
            file=str(input("Qual o nome do ficheiro? "))+".xml"
            parse_blast(file)
        elif ans=="14":
            ans = False
        else:
            print("\nInvalid")



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
    
