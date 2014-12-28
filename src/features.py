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


def features_CDS(record):
    featcds = [ ] 
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if record.features[i].type == "CDS": 
            featcds.append(i)
            if "locus_tag" in my_cds.qualifiers:
                print(my_cds.qualifiers["locus_tag"])
            else:
                print("Nao contem locus_tag!")
                
            if "product" in my_cds.qualifiers:
                print(my_cds.qualifiers["product"])
            else:
                print("Nao contem produtos!")
                
            if "note" in my_cds.qualifiers:
                print(my_cds.qualifiers["note"])
            else:
                print("Nao contem nota!")
    for k in featcds: 
        print (record.features[k].location)
        
def locus_tag(record):
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if record.features[i].type == "CDS": 
            featcds.append(i)
            if "locus_tag" in my_cds.qualifiers:
                print(my_cds.qualifiers["locus_tag"])
            else:
                print("Nao contem locus_tag!")
    

def gene_ID(record):
    for i in range(len(record.features)):
        my_gene = record.features[i]
        if record.features[i].type == "gene": 
            if "db_xref" in my_gene.qualifiers:
                print(my_gene.qualifiers["db_xref"])
 


#x=record.features[1]
#print(x.qualifiers)

#for feat in record.features: 
 #   print (str(feat))


def menu(record):
    ans=True
    while ans:
        print("""
    1.locus_tag
    2.CDS_product
    3.CDS_note
    4.geneID
    5.CDS_translation
    6.CDS_EC_number
    7.Sair
    """)
        ans=input("Qual a opcao? ")
        if ans=="1":
            features_CDS(record)
        #elif ans=="2":
            ######
        #elif ans=="3":
            #####
        elif ans=="4":
            gene_ID(record)
        #elif ans=="5":
            #####
        #elif ans=="6":
            #####
        elif ans=="7":
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
    



