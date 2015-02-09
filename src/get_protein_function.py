# -*- coding: utf-8 -*-

import os.path
import pandas #creating data frame
import numpy as np #creating array for data frame



#
#
#lista=[]        
#for name in range(len(GI)):
#    lista.append([])
#    lista[name].append(GI[name])
#    fic=open("../res/blast_with_note/match/allhits/funcao_all_hits/"+GI[name]+".txt").readlines()
#    #data = fic.replace("\n", " ")
#    for i in range(len(fic)):
#        data = fic[i].replace("\n", " ")
#        first = 'Gi: '
#        last = '  Possivel'
#        try:
#            start = data.rindex( first ) + len( first )
#            end = data.rindex(last, start)
#            novo= data[start:end]
#            lista[name].append(novo)
#        except:
#            pass       
#
#        first2 = 'função:   '
#        last2 = '{ECO:'
#        try:
#            start2 = data.rindex( first2 ) + len( first2 )
#            end2 = data.rindex(last2, start2)
#            novo2= data[start2:end2]
#            lista[name].append(novo2)
#        except:
#            pass       


#
#def tabela(l):
#    nome=[]
#    #dado=[]
#    for i in range(len(lista)):
#        #dado.append(lista[i][0])       
#        nome.append(lista[i])
#    data=np.array(nome)
#    df=pandas.DataFrame(data)
#    df.to_csv("../res/excel/funcao2", sep=',')
#
#print(tabela(lista))
def func():
    GI=[]
    for file in os.listdir("../res/blast_with_note/match/allhits/funcao_all_hits/"):
        if file.endswith(".txt"):
            limpo = file.replace(".txt", "")
            GI.append(limpo)
            
    lista=[]        
    for name in range(len(GI)):
        lista.append([])
        lista[name].append(GI[name])
        fic=open("../res/blast_with_note/match/allhits/funcao_all_hits/"+GI[name]+".txt").read()
        lista[name].append(fic)
    return lista
def tabela():
    lista=func()
    nome=[]
    #dado=[]
    for i in range(len(lista)):
        #dado.append(lista[i][0])       
        nome.append(lista[i])
    data=np.array(nome)
    df=pandas.DataFrame(data)
    df.to_csv("../res/excel/funcao2", sep=',')

tabela()