
 
from Bio.SwissProt import KeyWList

handle = open("uniprot.txt")
records = KeyWList.parse(handle)
x=[]
for record in records:
    x.append(record['ID'])
    
print(x)
handle.close()
