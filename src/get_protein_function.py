import os.path


GI=[]
for file in os.listdir("../res/blast_with_note/match/allhits/funcao_all_hits/"):
    if file.endswith(".txt"):
        limpo = file.replace(".txt", "")
        GI.append(limpo)

lista=[]        
for name in range(len(GI)):
    lista.append([])
    lista[name].append(GI[name])
    fic=open("../res/blast_with_note/match/allhits/funcao_all_hits/"+GI[name]+".txt").readlines()
    #data = fic.replace("\n", " ")
    for i in range(len(fic)):
        data = fic[i].replace("\n", " ")
        first = 'Gi: '
        last = '  Possivel'
        try:
            start = data.rindex( first ) + len( first )
            end = data.rindex(last, start)
            novo= data[start:end]
            lista[name].append(novo)
        except:
            pass       

    




