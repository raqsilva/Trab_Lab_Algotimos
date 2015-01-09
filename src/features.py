from Bio import SeqIO 
from Bio import Entrez
from Bio.Blast import NCBIWWW
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
    src = "D:\\Documentos\\GitHub\\Trab_Lab_Algotimos\\src\\"+filename #source folder
    dst = "D:\\Documentos\\GitHub\\Trab_Lab_Algotimos\\res"            #destination folder
    shutil.move(src, dst)
    record = SeqIO.read("../res/"+filename, "genbank") 
    return record


#Gives the sequence location in genome
def location(record):
    featcds = [ ] 
    for i in range(len(record.features)):
        if record.features[i].type == "CDS": 
            featcds.append(i)
    for k in featcds: 
        print (record.features[k].location)
        
        
#Notes from the sequence
def note(record):
     for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
            if "note" in my_cds.qualifiers:
                print(my_cds.qualifiers["note"])
            else:
                print("Nao contem nota!")
    
    
def locus_tag(record):
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS": 
            if "locus_tag" in my_cds.qualifiers:
                print(my_cds.qualifiers["locus_tag"])
            else:
                print("Nao contem locus_tag!")
    
    
#Products from the sequence, name of the proteins
def product(record):
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
            if "product" in my_cds.qualifiers:
                print(my_cds.qualifiers["product"])
            else:
                print("Nao contem produtos!")
       
         
#Gene ID and GI number
def gene_ID_GI(record):
    for i in range(len(record.features)):
        my_gene = record.features[i]
        if my_gene.type == "CDS": 
            #there are usually several db_xref entries
            if "db_xref" in my_gene.qualifiers:
                print(my_gene.qualifiers["db_xref"])
    
  
#protein EC number, identification
def EC_number(record):
     for i in range(len(record.features)):
        my_gene = record.features[i]
        if my_gene.type == "CDS": 
            if "EC_number" in my_gene.qualifiers:
                print(my_gene.qualifiers["EC_number"])
            else:
                print("Nao contem EC_number!")
                

#Translated sequence, protein sequence
def translation(record,gene):
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
            if my_cds.qualifiers["locus_tag"][0]==gene:
                if "translation" in my_cds.qualifiers:
                    return (my_cds.qualifiers["translation"])
                else:
                    return "Nao ha traducao"
            else:
                return "Nao existe este gene"
  

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
     
     
#Running Blast and saving info into a file
#GI_number - gene identification number
def blast(GI_numb,filename):
    result_handle = NCBIWWW.qblast("blastp", "swissprot", GI_numb)
    save_file = open(filename, "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()
    #moving the file to another directory
    src = "D:\\Documentos\\GitHub\\Trab_Lab_Algotimos\\src\\"+filename #source folder
    dst = "D:\\Documentos\\GitHub\\Trab_Lab_Algotimos\\res"            #destination folder
    shutil.move(src, dst)
    
    
#Parsing Blast files
def parse_blast(filename):
    result_handle = open("D:\\Documentos\\GitHub\\Trab_Lab_Algotimos\\res"+filename)



    
#x=record.features[1]
#print(x.qualifiers)

#for feat in record.features: 
 #   print (str(feat))


def menu(record):
    ans=True
    while ans:
        print("""
    1.locus_tag
    2.product
    3.note
    4.GI number, geneID
    5.translation
    6.EC_number
    7.location
    8.Artigos relacionados com um gene
    9.Correr o blast da proteina
    10.Parsing blast
    11.Sair
    """)
        ans=input("Qual a opcao? ")
        if ans=="1":
            locus_tag(record)
        elif ans=="2":
            product(record)
        elif ans=="3":
            note(record)
        elif ans=="4":
            gene_ID_GI(record)
        elif ans=="5":
            gene=str(input("Gene locus_tag: "))
            print(translation(record,gene))
        elif ans=="6":
            EC_number(record)
        elif ans=="7":
            location(record)
        elif ans=="8":
            gene=str(input("Gene name: "))
            print(DB_pubmed(gene))
        elif ans=="9":
            file=str(input("Qual o nome a colocar no ficheiro? "))+".xml"
            while os.path.isfile("D:\\Documentos\\GitHub\\Trab_Lab_Algotimos\\res\\"+file):
                 file=str(input("Qual o nome a colocar no ficheiro? "))+".xml"
            else:False
            GI=str(input("Qual o GI number da sequencia? "))
            blast(GI,file)
        elif ans=="10":
            file=str(input("Qual o nome do ficheiro? "))+".xml"
            parse_blast(file)
        elif ans=="11":
            ans = False
        else:
            print("\nInvalido")



def create_file():
    start=str(input(print("Insira o inicio da sua zona do genoma:")))
    stop=str(input(print("Insira o fim da sua zona do genoma:")))
    filename=str(input(print("Insira nome do ficheiro: "))+".gb")#filename plus extension being genbank
    #checking if name already exists    
    while os.path.isfile("D:\\Documentos\\GitHub\\Trab_Lab_Algotimos\\res\\"+filename):
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
    



