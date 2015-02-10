# -*- coding: utf-8 -*-
"""
Created on Tue Feb 10 14:18:05 2015

@author: Danielbraga
"""


from Bio import SeqIO #reading gb file 
from Bio.Blast import NCBIWWW,NCBIXML #fetching/parsing blast
import shutil#moving files
import os.path#cheking files in path
import urllib #getting info from site
from Uniprot_Parser import * #parsing uniprot text file



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
    for i in range(len(gi)):
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
        
#go to blast results and analise the hits, if e-value > 0.05 we don't consider this result        
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


#return all hits from blast and return one list with our protein and all aceptable hits from blast
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

#read allhits.txt, go to uniprot by hits and save txt by protein hit
def uniprotallhits():
    handle = open("../res/blast_without_note/match/allhits/allhits.txt").readlines()
    for n in range(len (handle)):       
            x=handle[n].split()
            for k in range(len(x)-1,len(x)):
                e=x[k]
                limpo=e.replace("']","")
                gi=(limpo[1:])
                #to organize results we create directory with gi of our proteins and save all txt from protein in the respectible gi
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


# go to foulder with all uniprot information from hits and and create a txt file with gi of our protein and all functions of hits
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
    1.Blast gi without note  
    2.Blast analiser without note
    3.Return all hits from blast without notes
    4.Go to uniprot and download information for all hits without notes
    5.Get all information for all hits without notes
    6.Exit
    """)
        ans=input("Choose an option? ")
        if ans=="1":   
            blastnote()
        elif ans=="2":
            blastanaliser()      
        elif ans=="3":
           file = open("../res/blast_without_note/match/allhits/allhits.txt",'w')
           lista=allhits()         
           for i in range(len(lista)):
               file.write("%s\n" % lista[i])
           file.close()
        elif ans=="4":
            uniprotallhits()
            
        elif ans=="5":    
            allfunctions()
        elif ans=="6":
            ans = False
        else:
            print("\nInvalid")

#main
if __name__ == "__main__":
    record = SeqIO.read("../res/"+'sequencia.gb', "genbank") 
    menu(record)