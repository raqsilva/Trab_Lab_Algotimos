handle = open("teste.txt").readlines()
text=str(handle)
string='SUBCELLULAR LOCATION:'
lista=[]
#for i in range (len(text)):
y=text.split(string)
m=y[1]
lista.append(m[0:10])
    
print(lista)