
 
#from Bio.SwissProt import KeyWList
from Bio.ExPASy import Enzyme


handle = open("uniprot.txt")
records = Enzyme.parse(handle)
x=[]
for record in records:
    x.append(record['CC'])
handle.close()
print(x)
#y=str(x[0])
##m=y.split("         ")
##print("nome",m[0],"revisao",m[1].replace)
#print (y)
