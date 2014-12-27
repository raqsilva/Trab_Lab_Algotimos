from Bio import SeqIO 
from Bio import Entrez
import shutil
import os.path

#[0:246000]

#fetching [0:246000] genome from Neisseria gonorrhoeae FA 1090 chromosome ncbi
#VERSION: NC_002946.2   GI:59800473

def get_genome_zone(start,stop,filename):
    #start=1
    #stop=246000
    #filename=filename +".gb"
    Entrez.email = "pg27668@alunos.uminho.pt"
    handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", seq_start="1", seq_stop="246000", id="59800473")
    file=open("sequencia.gb","w")#creating a GenBank file
    file.write(handle.read())
    file.close()
    handle.close()
    #moving the file to another directory
    src = "D:\\Documentos\\GitHub\\Trab_Lab_Algotimos\\src\\sequencia.gb"
    dst = "D:\\Documentos\\GitHub\\Trab_Lab_Algotimos\\res"
    shutil.move(src, dst)
    record = SeqIO.read("../res/sequencia.gb", "genbank") 
    return record


def features
#analysing the features
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


x=record.features[1]
print(x.qualifiers)


def teste():
    start=input(print("Insira o inicio da sua zona do genoma:"))
    stop=input(print("Insira o fim da sua zona do genoma:"))
    filename=input(print("Insira nome do ficheiro: "))
    while os.path.isfile("D:\\Documentos\\GitHub\\Trab_Lab_Algotimos\\src\\"+filename):
        filename=input(print("Insira outro nome: "))
        False
    record=get_genome_zone(start,stop,filename)
    


 
#if __name__ == "__main__":
#    teste()
    

#for feat in record.features: 
 #   print (str(feat))


