from Bio import SeqIO 
from Bio import Entrez
import shutil
import os.path

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
    src = "D:\\Documentos\\GitHub\\Trab_Lab_Algotimos\\src\\"+filename
    dst = "D:\\Documentos\\GitHub\\Trab_Lab_Algotimos\\res"
    shutil.move(src, dst)
    record = SeqIO.read("../res/"+filename, "genbank") 
    return record


def location(record):
    featcds = [ ] 
    for i in range(len(record.features)):
        if record.features[i].type == "CDS": 
            featcds.append(i)
    for k in featcds: 
        print (record.features[k].location)
        
        
def note(record):
     for i in range(len(record.features)):
        my_cds = record.features[i]
        if record.features[i].type == "CDS":
            if "note" in my_cds.qualifiers:
                print(my_cds.qualifiers["note"])
            else:
                print("Nao contem nota!")
    
    
def locus_tag(record):
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if record.features[i].type == "CDS": 
            if "locus_tag" in my_cds.qualifiers:
                print(my_cds.qualifiers["locus_tag"])
            else:
                print("Nao contem locus_tag!")
    

def product(record):
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if record.features[i].type == "CDS":
            if "product" in my_cds.qualifiers:
                print(my_cds.qualifiers["product"])
            else:
                print("Nao contem produtos!")
                

def gene_ID(record):
    for i in range(len(record.features)):
        my_gene = record.features[i]
        if record.features[i].type == "gene": 
            if "db_xref" in my_gene.qualifiers:
                print(my_gene.qualifiers["db_xref"])


def translation(record,gene):
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if record.features[i].type == "CDS":
            if my_cds.qualifiers["locus_tag"][0]==gene:
                if "translation" in my_cds.qualifiers:
                    return (my_cds.qualifiers["translation"])
                else:
                    return "Nao ha traducao"
            else:
                return "Nao existe este gene"
    

def EC_number(record):
     for i in range(len(record.features)):
        my_gene = record.features[i]
        if record.features[i].type == "CDS": 
            if "EC_number" in my_gene.qualifiers:
                print(my_gene.qualifiers["EC_number"])
            else:
                print("Nao contem EC_number!")
                


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
    4.geneID
    5.translation
    6.EC_number
    7.location
    8.Sair
    """)
        ans=input("Qual a opcao? ")
        if ans=="1":
            locus_tag(record)
        elif ans=="2":
            product(record)
        elif ans=="3":
            note(record)
        elif ans=="4":
            gene_ID(record)
        elif ans=="5":
            gene=str(input("Gene: "))
            print(translation(record,gene))
        elif ans=="6":
            EC_number(record)
        elif ans=="7":
            location(record)
        elif ans=="8":
            ans = False
        else:
            print("\nInvalido")


def create_file():
    start=str(input(print("Insira o inicio da sua zona do genoma:")))
    stop=str(input(print("Insira o fim da sua zona do genoma:")))
    filename=str(input(print("Insira nome do ficheiro: "))+".gb")
    while os.path.isfile("D:\\Documentos\\GitHub\\Trab_Lab_Algotimos\\res\\"+filename):
        filename=str(input(print("Insira outro nome: "))+".gb")
    else:False
    record=get_genome_zone(start,stop,filename)
    return record
    


#main
if __name__ == "__main__":
    res=str(input(print("Ja tem ficheiro S/N? "))).upper()
    if res=="N":
        record=create_file()
    else:
        filename=str(input("Nome do ficheiro? "))+".gb"
        record = SeqIO.read("../res/"+filename, "genbank") 
    menu(record)
    



