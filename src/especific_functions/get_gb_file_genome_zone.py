# -*- coding: utf-8 -*-
from Bio import SeqIO #reading gb file 
from Bio import Entrez #fetching genbank file
import shutil#moving files
import os.path#cheking files in path



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
    



