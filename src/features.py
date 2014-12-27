from Bio import SeqIO 
#[0:246000]

record = SeqIO.read("../res/sequence.gb", "genbank") 

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


 
#if __name__ == "__main__":
#    teste()
    

#for feat in record.features: 
 #   print (str(feat))

"""
def index_genbank_features(gb_record, feature_type, qualifier) :
    answer = dict()
    for (index, feature) in enumerate(gb_record.features) :
        if feature.type==feature_type :
            if qualifier in feature.qualifiers :
                #There should only be one locus_tag per feature, but there
                #are usually several db_xref entries
                for value in feature.qualifiers[qualifier] :
                    if value in answer :
                        print "WARNING - Duplicate key %s for %s features %i and %i" \
                           % (value, feature_type, answer[value], index)
                    else :
                        answer[value] = index
    return answer
"""

