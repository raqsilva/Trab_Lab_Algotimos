# -*- coding: utf-8 -*-e
import pandas #creating data frame
import numpy as np #creating array for data frame
from Uniprot_Info import *

def tabela(lista,locus):
    headers=['locus','gene name','location','GI GeneID','EC','proteinID','Uniprot_ID','product','note']
    dado=[]    
    for i in range(len(locus)):
        dado.append(locus[i])
    data=np.array(lista)
    df=pandas.DataFrame(data, dado, headers)
    df.to_csv("../res/excel/teste", sep='\t')
    

def tabela_pseudogenes(lista,locus):
    headers=['locus','location','gene_ID']
    dado=[]    
    for i in range(len(locus)):
        dado.append(locus[i])
    data=np.array(lista)
    df=pandas.DataFrame(data, dado, headers)
    df.to_csv("../res/excel/teste_pseudogenes", sep='\t')
    

def tabela_tRNA(lista,locus):
    headers=['locus','gene_ID','product','location']
    dado=[]    
    for i in range(len(locus)):
        dado.append(locus[i])
    data=np.array(lista)
    df=pandas.DataFrame(data, dado, headers)
    df.to_csv("../res/excel/teste_tRNA", sep='\t')
    

#csv file
def tabela_CDD(GI):
    dado=[]
    lista=[]
    for i in range(len(GI)):
        dado.append(GI[i][0])
        lista.append(GI[i])
    data=np.array(lista)
    df=pandas.DataFrame(data, dado)
    df.to_csv("../res/excel/CDD", sep='\t')
        

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

