from Bio import SeqIO 
 
record = SeqIO.read("../res/sequence.gb", "genbank") 

sub_record = record[0:246000]

featcds = [ ] 
for i in range(len(sub_record.features)):
    my_cds = sub_record.features[i]
    if sub_record.features[i].type == "CDS": 
        featcds.append(i)
        print(my_cds.qualifiers["locus_tag"])
        print(my_cds.qualifiers["product"])
        if "note" in my_cds.qualifiers:
            print(my_cds.qualifiers["note"])
        else:
            print("Nao contem nota!")
for k in featcds: 
    print (sub_record.features[k].location)


x=sub_record.features[1]
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

