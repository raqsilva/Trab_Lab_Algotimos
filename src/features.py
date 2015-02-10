# -*- coding: utf-8 -*-
from Bio import SeqIO #reading gb file 
from Bio import Entrez #fetching genbank file
from Bio.Blast import NCBIWWW,NCBIXML #fetching/parsing blast
import shutil#moving files
import os.path#cheking files in path
import urllib #getting info from site
from Uniprot_Parser import * #parsing uniprot text file
from Bio.SeqIO import UniprotIO #parsing uniprot xml file
import pandas #creating data frame
import numpy as np #creating array for data frame


#[0:246000] my zone

#fetching genome zone from Neisseria gonorrhoeae FA 1090 chromosome ncbi
#VERSION: NC_002946.2   GI:59800473


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


#Get all the info from one gene
def info(record,locus):
    lista=[]
    for i in range(len(locus)):
        lista.append([])
        lista[i].append(locus[i])
        lista[i].append(genes_names(record,locus[i]))
        lista[i].append(location(record,locus[i]))
        lista[i].append(gene_ID_GI(record,locus[i]))
        lista[i].append(EC_number(record,locus[i]))
        lista[i].append(protein_ID(record,locus[i]))
        lista[i].append(ID[i])
        lista[i].append(product(record,locus[i]))
        lista[i].append(note(record,locus[i]))
    return lista
    
    
def aceder(lista,nr):
    return lista[nr-1]


def tabela(lista,locus):
    headers=['locus','gene name','location','GI GeneID','EC','proteinID','Uniprot_ID','product','note']
    dado=[]    
    for i in range(len(locus)):
        dado.append(locus[i])
    data=np.array(lista)
    df=pandas.DataFrame(data, dado, headers)
    df.to_csv("../res/excel/teste", sep='\t')
    
    
#return protein ID
def protein_ID(record,locus_tag):
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
            if my_cds.qualifiers["locus_tag"][0]==str(locus_tag):
                if "protein_id" in my_cds.qualifiers:
                    return my_cds.qualifiers["protein_id"][0]
                else:
                    return "Nao contem protein_id"


#Gives the sequence location in genome
def location(record,locus_tag):
    featcds = [ ]
    for i in range(len(record.features)):
        if record.features[i].type == "CDS":
            if record.features[i].qualifiers["locus_tag"][0]==str(locus_tag): 
                featcds.append(i)
    for k in featcds: 
        return record.features[k].location
                
        
#Notes from the sequence
def note(record,locus_tag):
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
            if my_cds.qualifiers["locus_tag"][0]==str(locus_tag):  
                if "note" in my_cds.qualifiers:
                    return my_cds.qualifiers["note"][0]
                else:
                    return "Nao contem nota!"


#proteins without notes
def without_note(record):
    ID=[]
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
            if "note" not in my_cds.qualifiers:
                ID.append(my_cds.qualifiers["protein_id"][0])
    return ID


#locus tag, genes
def locus_tag(record):
    locus=[]
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS": 
            locus.append(my_cds.qualifiers["locus_tag"][0])
    return locus
    
    
#Products from the sequence, name of the proteins
def product(record,locus_tag):
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
            if my_cds.qualifiers["locus_tag"][0]==str(locus_tag):  
                if "product" in my_cds.qualifiers:
                    return my_cds.qualifiers["product"][0]
                else:
                    return "Nao contem produtos!"


#get hypothetical proteins (proteinID) from my zone
def hypoth_proteins(record):
    hypo=[]
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
            if my_cds.qualifiers["product"]==["hypothetical protein"]:
                hypo.append(my_cds.qualifiers["protein_id"][0])
    return hypo
    
       
#YP_009115478.1  
#Gene ID and GI number
def gene_ID_GI(record,locus_tag):
    for i in range(len(record.features)):
        my_gene = record.features[i]
        if my_gene.type == "CDS":             
        #there are usually several db_xref entries in CDS, but one in gene qualifier
            if my_gene.qualifiers["locus_tag"][0]==str(locus_tag):    
                if "db_xref" in my_gene.qualifiers:
                    x=my_gene.qualifiers["db_xref"]
                    return x[0],x[1]


def GI_number(record,locus_tag):
    GI=[]
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
            for j in range(len(locus_tag)):
                if my_cds.qualifiers["locus_tag"][0]==str(locus_tag[j]):
                    if "db_xref" in my_cds.qualifiers:
                        x=my_cds.qualifiers["db_xref"]
                        GI.append(x[0])
    return GI
 

#fetching CDD
def CDD(GI):
    CDD=[]
    for i in range(len(GI)):
        CDD.append([])
        CDD[i].append(GI[i])
        handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id=GI[i])
        record=SeqIO.read(handle, "gb")
        handle.close()
        for j in range(len(record.features)):
            my_cds = record.features[j]
            if my_cds.type == "Region" or my_cds.type == "Site":
                if "db_xref" in my_cds.qualifiers:
                    CDD[i].append(my_cds.qualifiers["db_xref"][0])
    return CDD    


GI=[['GI:59800474', 'CDD:256535', 'CDD:273037', 'CDD:99707', 'CDD:99707', 'CDD:99707', 'CDD:99707', 'CDD:99707', 'CDD:119330', 'CDD:119330'], ['GI:59800475', 'CDD:235541', 'CDD:238082', 'CDD:238082', 'CDD:238082', 'CDD:238082', 'CDD:238082'], ['GI:59800476', 'CDD:235469', 'CDD:257492', 'CDD:251338', 'CDD:197262', 'CDD:197262', 'CDD:197262', 'CDD:197262', 'CDD:197265', 'CDD:197265', 'CDD:197265', 'CDD:197265'], ['GI:59800477'], ['GI:59800478', 'CDD:234743', 'CDD:173906', 'CDD:173906', 'CDD:257918', 'CDD:275460', 'CDD:173912', 'CDD:173912', 'CDD:153412', 'CDD:153412'], ['GI:745998703'], ['GI:59800483', 'CDD:225495', 'CDD:238045', 'CDD:238045', 'CDD:238045', 'CDD:238045'], ['GI:59800485', 'CDD:252199'], ['GI:59800486', 'CDD:238190', 'CDD:238190', 'CDD:238190', 'CDD:238190'], ['GI:59800487', 'CDD:224867'], ['GI:59800488', 'CDD:225316', 'CDD:100107', 'CDD:100107'], ['GI:59800489', 'CDD:238786', 'CDD:238786'], ['GI:59800490', 'CDD:227114', 'CDD:238657', 'CDD:238657', 'CDD:238657'], ['GI:59800491', 'CDD:223734'], ['GI:59800492', 'CDD:223687', 'CDD:238560', 'CDD:238560'], ['GI:59800494', 'CDD:225117', 'CDD:249644'], ['GI:59800495', 'CDD:224422', 'CDD:249771'], ['GI:59800496', 'CDD:273856', 'CDD:239770', 'CDD:239770', 'CDD:239770', 'CDD:239770', 'CDD:173926', 'CDD:173926'], ['GI:59800497', 'CDD:264474'], ['GI:59800498', 'CDD:206754', 'CDD:206754'], ['GI:59800499', 'CDD:275600', 'CDD:198424'], ['GI:59800500', 'CDD:223532', 'CDD:173926', 'CDD:173926', 'CDD:276195'], ['GI:59800501', 'CDD:224135', 'CDD:250152'], ['GI:59800502', 'CDD:255039'], ['GI:59800503', 'CDD:238477', 'CDD:238477', 'CDD:238477', 'CDD:238477', 'CDD:238477'], ['GI:59800504', 'CDD:227307', 'CDD:145844', 'CDD:271179', 'CDD:271179'], ['GI:59800505', 'CDD:235470', 'CDD:238965', 'CDD:238965', 'CDD:132916', 'CDD:132916', 'CDD:132916', 'CDD:132916', 'CDD:251527'], ['GI:59800506', 'CDD:237673', 'CDD:250224', 'CDD:100105', 'CDD:100105', 'CDD:110897'], ['GI:59800508', 'CDD:234607', 'CDD:99735', 'CDD:99735', 'CDD:99735', 'CDD:99735'], ['GI:59800509', 'CDD:99838', 'CDD:99838', 'CDD:99838', 'CDD:99838', 'CDD:99838'], ['GI:59800510', 'CDD:253677', 'CDD:276194', 'CDD:100107', 'CDD:100107'], ['GI:59800511', 'CDD:236307', 'CDD:249744', 'CDD:190425', 'CDD:214878'], ['GI:59800512', 'CDD:235777', 'CDD:133459', 'CDD:133459', 'CDD:133459'], ['GI:59800513'], ['GI:59800514', 'CDD:234666'], ['GI:59800515', 'CDD:235393', 'CDD:249744', 'CDD:190425', 'CDD:198164', 'CDD:249744', 'CDD:276202', 'CDD:238712', 'CDD:238712', 'CDD:238712', 'CDD:238712', 'CDD:238712'], ['GI:59800516'], ['GI:59800517', 'CDD:238509', 'CDD:238509'], ['GI:59800518', 'CDD:237139', 'CDD:198165', 'CDD:153215', 'CDD:153215', 'CDD:153215'], ['GI:59800520', 'CDD:225953', 'CDD:253258'], ['GI:59800521', 'CDD:183270', 'CDD:238994', 'CDD:238994'], ['GI:59800522', 'CDD:239309', 'CDD:239309'], ['GI:59800523', 'CDD:188209', 'CDD:197670', 'CDD:238042', 'CDD:238042', 'CDD:238042'], ['GI:59800524', 'CDD:131349'], ['GI:59800525', 'CDD:224129', 'CDD:133044', 'CDD:133044', 'CDD:133044'], ['GI:59800526', 'CDD:238266', 'CDD:238266', 'CDD:238266'], ['GI:59800527', 'CDD:236583', 'CDD:206687', 'CDD:206687', 'CDD:206687', 'CDD:206687', 'CDD:206687', 'CDD:206687', 'CDD:206687', 'CDD:206687', 'CDD:206687', 'CDD:133440'], ['GI:59800528'], ['GI:59800529', 'CDD:224748', 'CDD:260768', 'CDD:238141'], ['GI:59800530', 'CDD:235645', 'CDD:173902', 'CDD:173902', 'CDD:173902', 'CDD:173902', 'CDD:173902'], ['GI:651851637', 'CDD:276322'], ['GI:59800531', 'CDD:235536', 'CDD:185679', 'CDD:185679', 'CDD:214901'], ['GI:59800532', 'CDD:235588', 'CDD:275460', 'CDD:173912', 'CDD:173912', 'CDD:173912', 'CDD:173909', 'CDD:173909', 'CDD:173909', 'CDD:153414', 'CDD:153414', 'CDD:153414', 'CDD:253933'], ['GI:59800533', 'CDD:276322'], ['GI:59800534', 'CDD:234739'], ['GI:59800535', 'CDD:260116', 'CDD:260116', 'CDD:260116', 'CDD:260116'], ['GI:59800536', 'CDD:223620', 'CDD:119389', 'CDD:119389', 'CDD:257599'], ['GI:59800539', 'CDD:235554', 'CDD:213988', 'CDD:213988', 'CDD:213988', 'CDD:254393', 'CDD:258717', 'CDD:239931', 'CDD:239931', 'CDD:239931'], ['GI:59800540', 'CDD:238095', 'CDD:238095'], ['GI:59800541'], ['GI:59800542', 'CDD:223686', 'CDD:213993', 'CDD:213993'], ['GI:59800543', 'CDD:236517', 'CDD:99734', 'CDD:99734', 'CDD:99734', 'CDD:99734'], ['GI:59800544', 'CDD:224011', 'CDD:271430', 'CDD:187548', 'CDD:187548', 'CDD:187548', 'CDD:187548', 'CDD:187548'], ['GI:59800545', 'CDD:99740', 'CDD:99740', 'CDD:99740', 'CDD:99740'], ['GI:59800546', 'CDD:251270', 'CDD:100050', 'CDD:234265', 'CDD:100050', 'CDD:100050'], ['GI:59800547', 'CDD:275912', 'CDD:223515'], ['GI:59800548', 'CDD:223515', 'CDD:275912'], ['GI:59800549', 'CDD:225153', 'CDD:263758'], ['GI:59800550', 'CDD:232920', 'CDD:238611', 'CDD:238611', 'CDD:238611', 'CDD:224896'], ['GI:59800551', 'CDD:234774', 'CDD:251987'], ['GI:59800552', 'CDD:238889', 'CDD:238889', 'CDD:249761', 'CDD:238889'], ['GI:59800553', 'CDD:173954', 'CDD:173954', 'CDD:173954', 'CDD:173954'], ['GI:59800554', 'CDD:223775', 'CDD:238260', 'CDD:238260', 'CDD:238260', 'CDD:238260'], ['GI:59800555', 'CDD:256588', 'CDD:233905', 'CDD:198033', 'CDD:252276', 'CDD:249725'], ['GI:59800556', 'CDD:225709'], ['GI:59800557', 'CDD:260846'], ['GI:59800558', 'CDD:225707'], ['GI:59800559', 'CDD:227306'], ['GI:59800560', 'CDD:227342', 'CDD:250219', 'CDD:250215'], ['GI:59800561', 'CDD:223296', 'CDD:206665', 'CDD:206665', 'CDD:206665', 'CDD:206665', 'CDD:206665', 'CDD:206665', 'CDD:206665', 'CDD:206665', 'CDD:206665'], ['GI:59800562', 'CDD:225418', 'CDD:249527'], ['GI:59800563', 'CDD:224252'], ['GI:59800564', 'CDD:275607'], ['GI:59800565', 'CDD:236585', 'CDD:276486'], ['GI:59800566', 'CDD:153246', 'CDD:153246'], ['GI:161572979', 'CDD:235374', 'CDD:249858', 'CDD:251520', 'CDD:111646'], ['GI:59800568', 'CDD:224938', 'CDD:251093'], ['GI:59800570'], ['GI:59800571', 'CDD:225407', 'CDD:252049'], ['GI:59800572', 'CDD:223715', 'CDD:119399', 'CDD:119399', 'CDD:119399', 'CDD:238030', 'CDD:238030', 'CDD:238030', 'CDD:238030'], ['GI:59800573', 'CDD:183285', 'CDD:239900', 'CDD:239900', 'CDD:239900', 'CDD:255792'], ['GI:59800574', 'CDD:239510', 'CDD:239510', 'CDD:239510'], ['GI:59800576', 'CDD:213561', 'CDD:238312', 'CDD:238312'], ['GI:59800577', 'CDD:236794', 'CDD:239934', 'CDD:239934', 'CDD:239934', 'CDD:238005', 'CDD:238005', 'CDD:238005', 'CDD:238034', 'CDD:238034', 'CDD:238034'], ['GI:59800578', 'CDD:234761', 'CDD:250373', 'CDD:187535', 'CDD:268288'], ['GI:59800579', 'CDD:179791'], ['GI:59800580'], ['GI:59800581', 'CDD:225406'], ['GI:59800582'], ['GI:59800583', 'CDD:130902', 'CDD:271629', 'CDD:271712', 'CDD:99778', 'CDD:99778', 'CDD:99778', 'CDD:99778', 'CDD:271757', 'CDD:257762'], ['GI:59800584'], ['GI:59800585'], ['GI:745998704', 'CDD:239100', 'CDD:239100', 'CDD:226359', 'CDD:276299', 'CDD:213179', 'CDD:213179', 'CDD:213179'], ['GI:59800586', 'CDD:179791', 'CDD:273260'], ['GI:59800588', 'CDD:238953', 'CDD:238953'], ['GI:59800589', 'CDD:148918'], ['GI:59800590', 'CDD:223792', 'CDD:250470', 'CDD:238351'], ['GI:59800591', 'CDD:226282', 'CDD:223675', 'CDD:100105', 'CDD:100105'], ['GI:59800593', 'CDD:235417'], ['GI:59800594', 'CDD:224671', 'CDD:252025'], ['GI:59800595'], ['GI:59800596', 'CDD:233695', 'CDD:257693', 'CDD:238487', 'CDD:238487', 'CDD:238487', 'CDD:238487'], ['GI:59800597', 'CDD:182661', 'CDD:238013', 'CDD:238013', 'CDD:238013', 'CDD:238013', 'CDD:238013', 'CDD:197771'], ['GI:59800598', 'CDD:224159'], ['GI:59800599'], ['GI:59800600', 'CDD:223809', 'CDD:119392', 'CDD:119392'], ['GI:59800601', 'CDD:224671', 'CDD:252025'], ['GI:59800602', 'CDD:226910', 'CDD:239963', 'CDD:215020'], ['GI:59800603', 'CDD:223396'], ['GI:59800604', 'CDD:234612', 'CDD:270364', 'CDD:270364', 'CDD:270364'], ['GI:651851638', 'CDD:257686'], ['GI:59800607', 'CDD:238167', 'CDD:214692', 'CDD:238167', 'CDD:238167', 'CDD:238167', 'CDD:238034', 'CDD:238034', 'CDD:238034'], ['GI:59800609', 'CDD:129820', 'CDD:239200', 'CDD:239200', 'CDD:239200', 'CDD:239200', 'CDD:239200'], ['GI:651851625', 'CDD:179348'], ['GI:161572978', 'CDD:238294', 'CDD:238294', 'CDD:238294', 'CDD:238294'], ['GI:59800612', 'CDD:153246', 'CDD:153246'], ['GI:59800614', 'CDD:99825', 'CDD:99825', 'CDD:99825', 'CDD:99825', 'CDD:255760'], ['GI:59800615', 'CDD:237585', 'CDD:189007', 'CDD:189007', 'CDD:189007'], ['GI:59800616'], ['GI:59800617'], ['GI:59800619'], ['GI:59800620', 'CDD:271796'], ['GI:59800621'], ['GI:59800622'], ['GI:59800624', 'CDD:223874', 'CDD:275504', 'CDD:238347'], ['GI:59800625', 'CDD:224033', 'CDD:275585', 'CDD:119348', 'CDD:119348', 'CDD:119348'], ['GI:59800626', 'CDD:224046', 'CDD:213202'], ['GI:59800627', 'CDD:235418'], ['GI:59800628', 'CDD:234581'], ['GI:59800629', 'CDD:223877', 'CDD:250860', 'CDD:253112'], ['GI:59800630', 'CDD:234588'], ['GI:59800631', 'CDD:223715', 'CDD:100122', 'CDD:100122', 'CDD:119399', 'CDD:119399', 'CDD:119399', 'CDD:238030', 'CDD:238030', 'CDD:238030', 'CDD:238030'], ['GI:59800632', 'CDD:223816', 'CDD:238088', 'CDD:238088', 'CDD:238088', 'CDD:238088', 'CDD:238088', 'CDD:238225', 'CDD:238225'], ['GI:59800633', 'CDD:225844', 'CDD:252887', 'CDD:256667'], ['GI:59800634', 'CDD:242419'], ['GI:59800635', 'CDD:238310', 'CDD:238310', 'CDD:238310'], ['GI:59800636', 'CDD:275617'], ['GI:59800637', 'CDD:234818', 'CDD:275655'], ['GI:59800638', 'CDD:235091'], ['GI:59800639', 'CDD:238607', 'CDD:238607', 'CDD:238607', 'CDD:238607'], ['GI:59800640', 'CDD:212141', 'CDD:212141'], ['GI:59800641', 'CDD:223991', 'CDD:176195', 'CDD:176195', 'CDD:176195'], ['GI:59800642', 'CDD:212541', 'CDD:212541', 'CDD:212541'], ['GI:59800643', 'CDD:251482'], ['GI:59800644', 'CDD:258018', 'CDD:235615', 'CDD:254281', 'CDD:273336'], ['GI:59800645', 'CDD:237275', 'CDD:254281', 'CDD:216990'], ['GI:59800646', 'CDD:238213', 'CDD:238213', 'CDD:238213', 'CDD:238213'], ['GI:59800647', 'CDD:183226', 'CDD:276299', 'CDD:213179', 'CDD:213179', 'CDD:213179', 'CDD:213179', 'CDD:213179', 'CDD:213179', 'CDD:213179', 'CDD:254779'], ['GI:59800648', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394'], ['GI:59800649', 'CDD:224098', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394'], ['GI:59800650', 'CDD:276298', 'CDD:223737'], ['GI:59800651', 'CDD:223083'], ['GI:59800652', 'CDD:236490', 'CDD:198027', 'CDD:239906', 'CDD:239906', 'CDD:238548', 'CDD:238548', 'CDD:238548', 'CDD:238548', 'CDD:238548'], ['GI:59800653', 'CDD:235809', 'CDD:250534', 'CDD:249825', 'CDD:271901'], ['GI:59800654', 'CDD:226201'], ['GI:59800655', 'CDD:224719'], ['GI:59800656', 'CDD:276303', 'CDD:223620', 'CDD:119389', 'CDD:119389'], ['GI:59800657', 'CDD:153097', 'CDD:153097'], ['GI:59800658', 'CDD:178807'], ['GI:59800659', 'CDD:270377', 'CDD:270377', 'CDD:270377'], ['GI:59800660', 'CDD:223539', 'CDD:132997', 'CDD:132997'], ['GI:59800661', 'CDD:223735', 'CDD:100051', 'CDD:100051', 'CDD:100051'], ['GI:59800662', 'CDD:179105', 'CDD:206732', 'CDD:206732', 'CDD:206732', 'CDD:206732', 'CDD:206732', 'CDD:206732', 'CDD:206732', 'CDD:206732', 'CDD:206732', 'CDD:206732', 'CDD:239660', 'CDD:258632'], ['GI:59800663', 'CDD:234598'], ['GI:59800664', 'CDD:240082', 'CDD:240082', 'CDD:240082'], ['GI:59800665', 'CDD:240083', 'CDD:240083'], ['GI:59800666', 'CDD:153219', 'CDD:171871', 'CDD:153219', 'CDD:153219', 'CDD:153219'], ['GI:59800667', 'CDD:223357'], ['GI:59800668', 'CDD:226361', 'CDD:213226', 'CDD:213226', 'CDD:213226', 'CDD:213226', 'CDD:213226', 'CDD:213226', 'CDD:213226', 'CDD:213226'], ['GI:59800669', 'CDD:224099', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394'], ['GI:59800670', 'CDD:270261', 'CDD:270261'], ['GI:59800672', 'CDD:179143', 'CDD:176463', 'CDD:176463', 'CDD:176463'], ['GI:59800673', 'CDD:133021', 'CDD:182252', 'CDD:133021', 'CDD:133021'], ['GI:59800674', 'CDD:238285', 'CDD:238285', 'CDD:238285'], ['GI:59800675'], ['GI:59800676', 'CDD:238239', 'CDD:238239', 'CDD:238239', 'CDD:238239'], ['GI:59800677', 'CDD:240022', 'CDD:240022'], ['GI:59800678', 'CDD:191481', 'CDD:258590'], ['GI:59800679'], ['GI:59800680', 'CDD:257184'], ['GI:59800681'], ['GI:59800682'], ['GI:59800683', 'CDD:223246', 'CDD:266653'], ['GI:59800684', 'CDD:234673', 'CDD:163665', 'CDD:163665', 'CDD:163665'], ['GI:59800685', 'CDD:100004', 'CDD:100004', 'CDD:100004'], ['GI:59800686', 'CDD:111368'], ['GI:59800687', 'CDD:236346', 'CDD:100105', 'CDD:100105', 'CDD:253997'], ['GI:59800688', 'CDD:236137', 'CDD:238062', 'CDD:238062', 'CDD:238062', 'CDD:238062', 'CDD:145978', 'CDD:251739', 'CDD:257328', 'CDD:237994', 'CDD:237994', 'CDD:237994'], ['GI:59800689', 'CDD:225629', 'CDD:129010'], ['GI:59800690', 'CDD:225567', 'CDD:133475', 'CDD:133475', 'CDD:133475', 'CDD:133475'], ['GI:59800692', 'CDD:234814', 'CDD:238835', 'CDD:238835', 'CDD:238835'], ['GI:59800693', 'CDD:223358', 'CDD:249824', 'CDD:133453', 'CDD:133453'], ['GI:59800695', 'CDD:234808']]


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
        

#Getting pseudogenes locus_tag
def pseudogenes(record):
    locus=[]
    for i in range(len(record.features)):
        my_gene = record.features[i]
        if my_gene.type == "gene":
             if "pseudogene" in my_gene.qualifiers:
                locus.append(my_gene.qualifiers["locus_tag"][0])
    return locus


def location_pseudo(record,locus_tag):
    featcds = [ ]
    for i in range(len(record.features)):
        if record.features[i].type == "gene":
            if record.features[i].qualifiers["locus_tag"][0]==str(locus_tag): 
                featcds.append(i)
    for k in featcds: 
        return record.features[k].location


def pseudogeneID(record,locus_tag):
    for i in range(len(record.features)):
        my_gene = record.features[i]
        if my_gene.type == "gene":
            if my_gene.qualifiers["locus_tag"][0]==str(locus_tag): 
                if "pseudogene" in my_gene.qualifiers:
                    return my_gene.qualifiers["db_xref"][0]


def info_pseudogenes(record,locus):
    lista=[]
    for i in range(len(locus)):
        lista.append([])
        lista[i].append(locus[i])
        lista[i].append(location_pseudo(record,locus[i]))
        lista[i].append(pseudogeneID(record,locus[i]))
    return lista   
    

def tabela_pseudogenes(lista,locus):
    headers=['locus','location','gene_ID']
    dado=[]    
    for i in range(len(locus)):
        dado.append(locus[i])
    data=np.array(lista)
    df=pandas.DataFrame(data, dado, headers)
    df.to_csv("../res/excel/teste_pseudogenes", sep='\t')
    

#Getting genes names
def genes_names(record,locus_tag):
    for i in range(len(record.features)):
        my_gene = record.features[i]
        if my_gene.type == "CDS":
            if my_gene.qualifiers["locus_tag"][0]==str(locus_tag):  
                if "gene" in my_gene.qualifiers:
                    return my_gene.qualifiers["gene"][0]
                else:
                    return "Nao tem nome!"


#list of genes names
def list_genes_names(record,locus):
    genes=[]
    for i in range(len(locus)):
        if genes_names(record,locus[i])!="Nao tem nome!":
            genes.append(genes_names(record,locus[i]))
    return genes
    
  
#protein EC number, identification
def EC_number(record,locus_tag):
    for i in range(len(record.features)):
        my_gene = record.features[i]
        if my_gene.type == "CDS":
            if my_gene.qualifiers["locus_tag"][0]==str(locus_tag):  
                if "EC_number" in my_gene.qualifiers:
                    return my_gene.qualifiers["EC_number"][0]
                else:
                    return "Nao contem EC_number!"
      

#Translated sequence, protein sequence
def translation(record,protein_ID): 
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
             if my_cds.qualifiers["protein_id"][0]==str(protein_ID):
                    if "translation" in my_cds.qualifiers:
                        print (my_cds.qualifiers["translation"][0])
                    else:
                        print( "Nao ha traducao")
                    

#locus from tRNA             
def tRNA(record):
    #trna=[]
    locus=[]
    for i in range(len(record.features)):
        my_trna = record.features[i]
        if my_trna.type == "tRNA":
            #trna.append(my_trna.qualifiers["product"][0])
            locus.append(my_trna.qualifiers["locus_tag"][0])
    return locus


#gets geneID form tRNA
def gene_ID_tRNA(record,locus_tag):
    for i in range(len(record.features)):
        my_gene = record.features[i]
        if my_gene.type == "tRNA":             
            if my_gene.qualifiers["locus_tag"][0]==str(locus_tag):    
                x=my_gene.qualifiers["db_xref"]
                return x[0]
                    

#gets product from tRNA
def product_tRNA(record,locus_tag):
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "tRNA":
            if my_cds.qualifiers["locus_tag"][0]==str(locus_tag):  
                return my_cds.qualifiers["product"][0]
                    

#Gets location strand from tRNA
def location_tRNA(record,locus_tag):
    featcds = [ ]
    for i in range(len(record.features)):
        if record.features[i].type == "tRNA":
            if record.features[i].qualifiers["locus_tag"][0]==str(locus_tag): 
                featcds.append(i)
    for k in featcds: 
        return record.features[k].location


#Gets all the info from tRNA                
def info_tRNA(record,locus):
    lista=[]
    for i in range(len(locus)):
        lista.append([])
        lista[i].append(locus[i])
        lista[i].append(gene_ID_tRNA(record,locus[i]))
        lista[i].append(product_tRNA(record,locus[i]))
        lista[i].append(location_tRNA(record,locus[i]))
    return lista   



def tabela_tRNA(lista,locus):
    headers=['locus','gene_ID','product','location']
    dado=[]    
    for i in range(len(locus)):
        dado.append(locus[i])
    data=np.array(lista)
    df=pandas.DataFrame(data, dado, headers)
    df.to_csv("../res/excel/teste_tRNA", sep='\t')
    
              

def uniprot_ID(record):
    locus=locus_tag(record)
    proteins=[]
    IDs=[]
    for i in range(len(locus)):
        proteins.append(protein_ID(record,locus[i]))
    for j in range(len(proteins)):
        data = urllib.request.urlopen("http://www.uniprot.org/uniprot/?query="+proteins[j]+"&sort=score").read()
        if len(data.split())>=3298:
            x=str(data.split()[3298])
            IDs.append(x[6:12])
        else:
            IDs.append(proteins[j])
    return IDs


#['b\'id="Q5FAJ2"\'', 'b\'id="Q5FAJ1"\'', 'b\'id="Q5FAJ0"\'', 'b\'id="Q5FAI9"\'', 'b\'id="Q5FAJ3"\'', 'YP_009115477.1', 'b\'id="Q5FAK0"\'', 'b\'id="Q5FAJ8"\'', 'b\'id="Q5FAJ7"\'', 'b\'id="Q5FAJ6"\'', 'b\'id="Q5FAJ5"\'', 'b\'id="Q5FAJ4"\'', 'b\'id="Q5FAL3"\'', 'b\'id="Q5FAL2"\'', 'b\'id="Q5FAL1"\'', 'b\'id="Q5FAK9"\'', 'b\'id="Q5FAK8"\'', 'b\'id="Q5FAK7"\'', 'b\'id="Q5FAK6"\'', 'b\'id="Q5FAK5"\'', 'b\'id="Q5FAI8"\'', 'b\'id="Q5FAI7"\'', 'b\'id="Q5FAI6"\'', 'b\'id="Q5FAI5"\'', 'b\'id="Q5FAI4"\'', 'b\'id="Q5FAI3"\'', 'b\'id="Q5FAI2"\'', 'b\'id="Q5FAI1"\'', 'b\'id="Q5FAH9"\'', 'b\'id="Q5FAH8"\'', 'b\'id="Q5FAH7"\'', 'b\'id="Q5FAH6"\'', 'b\'id="Q5FAH5"\'', 'b\'id="Q5FAH4"\'', 'b\'id="Q5FAH3"\'', 'b\'id="Q5FAH2"\'', 'b\'id="Q5FAH1"\'', 'b\'id="Q5FAH0"\'', 'b\'id="Q5FAG9"\'', 'b\'id="Q5FAG7"\'', 'b\'id="Q5FAG6"\'', 'b\'id="Q5FAG5"\'', 'b\'id="Q5FAG4"\'', 'b\'id="Q5FAG3"\'', 'b\'id="Q5FAG2"\'', 'b\'id="Q5FAG1"\'', 'b\'id="Q5FAG0"\'', 'b\'id="Q5FAF9"\'', 'b\'id="Q5FAF8"\'', 'b\'id="Q5FAF7"\'', 'YP_008914846.1', 'b\'id="Q5FAF6"\'', 'b\'id="Q5FAF5"\'', 'b\'id="Q5FAF4"\'', 'b\'id="Q5FAF3"\'', 'b\'id="Q5FAF2"\'', 'b\'id="Q5FAF1"\'', 'b\'id="Q5FAE8"\'', 'b\'id="Q5FAE7"\'', 'b\'id="Q5FAE6"\'', 'b\'id="Q5FAE5"\'', 'b\'id="Q5FAE4"\'', 'b\'id="Q5FAE3"\'', 'b\'id="Q5FAE2"\'', 'b\'id="Q5FAE1"\'', 'b\'id="Q5FAE0"\'', 'b\'id="Q5FAD9"\'', 'b\'id="Q5FAD8"\'', 'b\'id="Q5FAD7"\'', 'b\'id="Q5FAD6"\'', 'b\'id="Q5FAD5"\'', 'b\'id="Q5FAD4"\'', 'b\'id="Q5FAD3"\'', 'b\'id="Q5FAD2"\'', 'b\'id="Q5FAD1"\'', 'b\'id="Q5FAD0"\'', 'b\'id="Q5FAC9"\'', 'b\'id="Q5FAC8"\'', 'b\'id="Q5FAC7"\'', 'b\'id="Q5FAC6"\'', 'b\'id="Q5FAC5"\'', 'b\'id="Q5FAC4"\'', 'b\'id="Q5FAC3"\'', 'b\'id="Q5FAC2"\'', 'b\'id="Q5FAC1"\'', 'b\'id="Q5FAC0"\'', 'b\'id="Q5FAB9"\'', 'b\'id="Q5FAB7"\'', 'b\'id="Q5FAB6"\'', 'b\'id="Q5FAB5"\'', 'b\'id="Q5FAB4"\'', 'b\'id="Q5FAB3"\'', 'b\'id="Q5FAB1"\'', 'b\'id="Q5FAB0"\'', 'b\'id="Q5FAA9"\'', 'b\'id="Q5FAA8"\'', 'b\'id="Q5FAA7"\'', 'b\'id="Q5FAA6"\'', 'b\'id="Q5FAA5"\'', 'b\'id="Q5FAA4"\'', 'b\'id="Q5FAA3"\'', 'b\'id="Q5FAA2"\'', 'YP_009115478.1', 'b\'id="Q5FAA1"\'', 'b\'id="Q5FA99"\'', 'b\'id="Q5FA98"\'', 'b\'id="Q5FA97"\'', 'b\'id="Q5FA96"\'', 'b\'id="Q5FA94"\'', 'b\'id="Q5FA93"\'', 'b\'id="Q5FA92"\'', 'b\'id="Q5FA91"\'', 'b\'id="Q5FA90"\'', 'b\'id="Q5FA89"\'', 'b\'id="Q5FA88"\'', 'b\'id="Q5FA87"\'', 'b\'id="Q5FA86"\'', 'b\'id="Q5FA85"\'', 'b\'id="Q5FA84"\'', 'b\'id="Q5FA83"\'', 'YP_008914847.1', 'b\'id="Q5FA80"\'', 'b\'id="Q5FA78"\'', 'YP_207322.2', 'b\'id="Q5FA76"\'', 'b\'id="Q5FA75"\'', 'b\'id="Q5FA73"\'', 'b\'id="Q5FA72"\'', 'b\'id="Q5FA71"\'', 'b\'id="Q5FA70"\'', 'b\'id="Q5FA68"\'', 'b\'id="Q5FA67"\'', 'b\'id="Q5FA66"\'', 'b\'id="Q5FA65"\'', 'b\'id="Q5FA63"\'', 'b\'id="Q5FA62"\'', 'b\'id="Q5FA61"\'', 'b\'id="Q5FA60"\'', 'b\'id="Q5FA59"\'', 'b\'id="Q5FA58"\'', 'b\'id="Q5FA57"\'', 'b\'id="Q5FA56"\'', 'b\'id="Q5FA55"\'', 'b\'id="Q5FA54"\'', 'b\'id="Q5FA53"\'', 'b\'id="Q5FA52"\'', 'b\'id="Q5FA51"\'', 'b\'id="Q5FA50"\'', 'b\'id="Q5FA49"\'', 'b\'id="Q5FA48"\'', 'b\'id="Q5FA47"\'', 'b\'id="Q5FA46"\'', 'b\'id="Q5FA45"\'', 'b\'id="Q5FA44"\'', 'b\'id="Q5FA43"\'', 'b\'id="Q5FA42"\'', 'b\'id="Q5FA41"\'', 'b\'id="Q5FA40"\'', 'b\'id="Q5FA39"\'', 'b\'id="Q5FA38"\'', 'b\'id="Q5FA37"\'', 'b\'id="Q5FA36"\'', 'b\'id="Q5FA35"\'', 'b\'id="Q5FA34"\'', 'b\'id="Q5FA33"\'', 'b\'id="Q5FA32"\'', 'b\'id="Q5FA31"\'', 'b\'id="Q5FA30"\'', 'b\'id="Q5FA29"\'', 'b\'id="Q5FA28"\'', 'b\'id="Q5FA27"\'', 'b\'id="Q5FA26"\'', 'b\'id="Q5FA25"\'', 'b\'id="Q5FA24"\'', 'b\'id="Q5FA23"\'', 'b\'id="Q5FA22"\'', 'b\'id="Q5FA21"\'', 'b\'id="Q5FA20"\'', 'b\'id="Q5FA19"\'', 'b\'id="Q5FA18"\'', 'b\'id="Q5FA17"\'', 'b\'id="Q5FA15"\'', 'b\'id="Q5FA14"\'', 'b\'id="Q5FA13"\'', 'b\'id="Q5FA12"\'', 'b\'id="Q5FA11"\'', 'b\'id="Q5FA10"\'', 'b\'id="Q5FA09"\'', 'b\'id="Q5FA08"\'', 'b\'id="Q5FA07"\'', 'b\'id="Q5FA06"\'', 'b\'id="Q5FA05"\'', 'b\'id="Q5FA04"\'', 'b\'id="Q5FA03"\'', 'b\'id="Q5FA02"\'', 'b\'id="Q5FA01"\'', 'b\'id="Q5FA00"\'', 'b\'id="Q5F9Z9"\'', 'b\'id="Q5F9Z8"\'', 'b\'id="Q5F9Z7"\'', 'b\'id="Q5F9Z5"\'', 'b\'id="Q5F9Z4"\'', 'b\'id="Q5F9Z2"\'']

ID=['Q5FAJ2', 'Q5FAJ1', 'Q5FAJ0', 'Q5FAI9', 'Q5FAJ3', 'YP_009115477.1', 'Q5FAK0', 'Q5FAJ8', 'Q5FAJ7', 'Q5FAJ6', 'Q5FAJ5', 'Q5FAJ4', 'Q5FAL3', 'Q5FAL2', 'Q5FAL1', 'Q5FAK9', 'Q5FAK8', 'Q5FAK7', 'Q5FAK6', 'Q5FAK5', 'Q5FAI8', 'Q5FAI7', 'Q5FAI6', 'Q5FAI5', 'Q5FAI4', 'Q5FAI3', 'Q5FAI2', 'Q5FAI1', 'Q5FAH9', 'Q5FAH8', 'Q5FAH7', 'Q5FAH6', 'Q5FAH5', 'Q5FAH4', 'Q5FAH3', 'Q5FAH2', 'Q5FAH1', 'Q5FAH0', 'Q5FAG9', 'Q5FAG7', 'Q5FAG6', 'Q5FAG5', 'Q5FAG4', 'Q5FAG3', 'Q5FAG2', 'Q5FAG1', 'Q5FAG0', 'Q5FAF9', 'Q5FAF8', 'Q5FAF7', 'YP_008914846.1', 'Q5FAF6', 'Q5FAF5', 'Q5FAF4', 'Q5FAF3', 'Q5FAF2', 'Q5FAF1', 'Q5FAE8', 'Q5FAE7', 'Q5FAE6', 'Q5FAE5', 'Q5FAE4', 'Q5FAE3', 'Q5FAE2', 'Q5FAE1', 'Q5FAE0', 'Q5FAD9', 'Q5FAD8', 'Q5FAD7', 'Q5FAD6', 'Q5FAD5', 'Q5FAD4', 'Q5FAD3', 'Q5FAD2', 'Q5FAD1', 'Q5FAD0', 'Q5FAC9', 'Q5FAC8', 'Q5FAC7', 'Q5FAC6', 'Q5FAC5', 'Q5FAC4', 'Q5FAC3', 'Q5FAC2', 'Q5FAC1', 'Q5FAC0', 'Q5FAB9', 'Q5FAB7', 'Q5FAB6', 'Q5FAB5', 'Q5FAB4', 'Q5FAB3', 'Q5FAB1', 'Q5FAB0', 'Q5FAA9', 'Q5FAA8', 'Q5FAA7', 'Q5FAA6', 'Q5FAA5', 'Q5FAA4', 'Q5FAA3', 'Q5FAA2', 'YP_009115478.1', 'Q5FAA1', 'Q5FA99', 'Q5FA98', 'Q5FA97', 'Q5FA96', 'Q5FA94', 'Q5FA93', 'Q5FA92', 'Q5FA91', 'Q5FA90', 'Q5FA89', 'Q5FA88', 'Q5FA87', 'Q5FA86', 'Q5FA85', 'Q5FA84', 'Q5FA83', 'YP_008914847.1', 'Q5FA80', 'Q5FA78', 'YP_207322.2', 'Q5FA76', 'Q5FA75', 'Q5FA73', 'Q5FA72', 'Q5FA71', 'Q5FA70', 'Q5FA68', 'Q5FA67', 'Q5FA66', 'Q5FA65', 'Q5FA63', 'Q5FA62', 'Q5FA61', 'Q5FA60', 'Q5FA59', 'Q5FA58', 'Q5FA57', 'Q5FA56', 'Q5FA55', 'Q5FA54', 'Q5FA53', 'Q5FA52', 'Q5FA51', 'Q5FA50', 'Q5FA49', 'Q5FA48', 'Q5FA47', 'Q5FA46', 'Q5FA45', 'Q5FA44', 'Q5FA43', 'Q5FA42', 'Q5FA41', 'Q5FA40', 'Q5FA39', 'Q5FA38', 'Q5FA37', 'Q5FA36', 'Q5FA35', 'Q5FA34', 'Q5FA33', 'Q5FA32', 'Q5FA31', 'Q5FA30', 'Q5FA29', 'Q5FA28', 'Q5FA27', 'Q5FA26', 'Q5FA25', 'Q5FA24', 'Q5FA23', 'Q5FA22', 'Q5FA21', 'Q5FA20', 'Q5FA19', 'Q5FA18', 'Q5FA17', 'Q5FA15', 'Q5FA14', 'Q5FA13', 'Q5FA12', 'Q5FA11', 'Q5FA10', 'Q5FA09', 'Q5FA08', 'Q5FA07', 'Q5FA06', 'Q5FA05', 'Q5FA04', 'Q5FA03', 'Q5FA02', 'Q5FA01', 'Q5FA00', 'Q5F9Z9', 'Q5F9Z8', 'Q5F9Z7', 'Q5F9Z5', 'Q5F9Z4', 'Q5F9Z2']


def info_uniprot():
    identifier=[]
    handle = open("../res/"+"uniprot.txt")
    records = parse(handle) # Uses the function 'parse' from the module. 
    for record in records:
        for i in range(len(ID)):
            identifier.append([])
            if record["AC"]==ID[i]+";":
                #identifier.append([])
                identifier[i].append(ID[i])
                identifier[i].append(record["ID"])
                identifier[i].append(record["DE"])
                identifier[i].append(record["CC"])
    handle.close()
    return identifier


def more_info_uniprot():
    handle = open("../res/"+"uniprot_XML.xml")
    records=UniprotIO.UniprotIterator(handle,return_raw_comments=True)
    refs=[]
    for record in records:
        for i in range(len(ID)):
            if record.id==ID[i]:
                refs.append([ID[i]]+record.dbxrefs)#GO´s
    handle.close()
    return refs
    

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


def GO():
    l=[]
    lista=sorting(tab())
    for i in range(len(lista)):
        for j in range(len(lista[i])):
            if "GO:GO:" in lista[i][j]:
                l.append(lista[i][0])
                l.append(lista[i][j])
    return l


#Searching articles from PubMed DB referring to my organism and a gene
def DB_pubmed(gene):
    """
    handle = Entrez.egquery(term="Neisseria gonorrhoeae")
    record = Entrez.read(handle)
    for row in record["eGQueryResult"]:
        if row["DbName"]=="pubmed":
            print(row["Count"])
    """
    handle = Entrez.esearch(db="pubmed", term="Neisseria gonorrhoeae[Orgn] AND "+gene+"[Gene]", retmax=11117)
    record = Entrez.read(handle)
    idlist = record["IdList"]
    return idlist


#745998704  
#Running Blast and saving info into a file
#GI_number - gene identification number
def blast(GI_numb,filename):
    result_handle = NCBIWWW.qblast("blastp", "swissprot", GI_numb)
    save_file = open(filename, "w")
    save_file.write(result_handle.read())
    save_file.close()
    result_handle.close()
    #moving the file to another directory
    path=os.getcwd()
    src = path+"/"+filename #source folder
    dst = "../res"#destination folder
    shutil.move(src, dst)
    
    
#Parsing Blast files
def parse_blast(filename):
    #E_VALUE_THRESH = 0.05
    result_handle = open("../res/"+filename)
    blast_record = NCBIXML.read(result_handle)
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            #if hsp.expect < E_VALUE_THRESH:
            print ("****Alignment****")
            print ('sequence:', alignment.title)
            print ('length:', alignment.length)
            print ('e value:', hsp.expect)
    result_handle.close()


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
    for i in range(6,len(gi)):
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

def gimatch():
    handle = open("../res/blast_without_note/match/matches.txt").readlines()    
    gimatch=[]
    gi=[]
    hit=[]
    s2='sp|'
    xml='['
    for i in handle:
        n = i[i.index(xml) + len(xml):] 
        g = n.split('.xml', 1)[0]  
        gi.append(g[1:])
    for J in handle:
        s3 = J[J.index(s2) + len(s2):] 
        sP = s3.split('|', 1)[0] 
        sp=sP[:6]    
        hit.append(sp)
    for k  in range(len(gi)):
        gimatch.append(gi[k]+' '+hit[k])    
            
            
    
    for n in range(len(gi)):
        site = urllib.request.urlopen("http://www.uniprot.org/uniprot/"+hit[n]+".txt")
        data = site.readlines()
        nome=gi[n]
        file = open("../res/blast_without_note/match/function/"+nome+'.txt',"wb") #open file in binary mode
        file.writelines(data)
        file.close()
        
        
#function to see possible function of proteins that didnt have note:
def getfunction():
    blast=[]    
    for file in os.listdir("../res/blast_without_note/match/function/"):
        if file.endswith(".txt"):
               blast.append(file)
    
    lista=[]
    for j in range(len(blast)):
        lista.append([])
        nome=blast[j]
        gi = nome.replace(".txt","")
        first = '-!- FUNCTION:'
        last = 'CC'
        file = open("../res/blast_without_note/match/function/"+nome).read()
        data = file.replace("\n", " ") 
        try:
            start = data.rindex( first ) + len( first )
            end = data.rindex(last, start)
            novo= data[start:end] 
            lista[j].append('Gi: '+gi+'  '+'Possivel função:  '+novo)
        except:
            pass   
    return lista
    
#return all hits from blast
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

def uniprotallhits():
    handle = open("../res/blast_without_note/match/allhits/allhits.txt").readlines()
    
    
    
    for n in range(len (handle)):       
            x=handle[n].split()
            for k in range(len(x)-1,len(x)):
                e=x[k]
                limpo=e.replace("']","")
                gi=(limpo[1:])
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


#GI numbers from genes with note
def GInumbers(record,locus_tag):
    GI=[]
    for i in range(len(record.features)):
        my_cds = record.features[i]
        if my_cds.type == "CDS":
            for j in range(len(locus_tag)):
                if my_cds.qualifiers["locus_tag"][0]==str(locus_tag[j]):
                    if "note" in my_cds.qualifiers:
                        if "db_xref" in my_cds.qualifiers:
                            x=my_cds.qualifiers["db_xref"]
                            GI.append(x[0])
    return GI

#blast gi with note    
def blastwithnote():
    locus=locus_tag(record)
    gi=GInumbers(record,locus)
    for i in range(37,len(gi)):
        gis=str(gi[i])
        GI_numb=gis[3:]
        result_handle = NCBIWWW.qblast("blastp","swissprot", GI_numb)
        save_file = open(GI_numb+'.xml', "w")
        save_file.write(result_handle.read())
        save_file.close()
        result_handle.close()
        #moving the file to another directory
        path=os.getcwd()
        src = path+"/"+GI_numb+'.xml' #source folder
        dst = "../res/blast_with_note/"#destination folder
        shutil.move(src, dst)
                              
#return all hits from blast
def allhitswithnote():
    lista=[]
    handle = open("../res/blast_with_note/match/matches.txt").readlines()    
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
 

lista_GO=[['Q5FAJ2', 'GO:GO:0003688', 'GO:GO:0005524', 'GO:GO:0005737',  'GO:GO:0006270', 'GO:GO:0006275'] ,['Q5FAJ1', 'GO:GO:0003677', 'GO:GO:0003887',  'GO:GO:0005737',  'GO:GO:0006260', 'GO:GO:0008408',  'GO:GO:0009360'],[ 'Q5FAJ0', 'GO:GO:0005524', 'GO:GO:0006799',  'GO:GO:0008976',  'GO:GO:0009358'],['Q5FAJ3', 'GO:GO:0002161', 'GO:GO:0004823',  'GO:GO:0005524', 'GO:GO:0005737', 'GO:GO:0006429'], ['Q5FAK0', 'GO:GO:0043565'],['Q5FAJ8', 'GO:GO:0009306', 'GO:GO:0015450',  'GO:GO:0016021'],['Q5FAJ7', 'GO:GO:0004807', 'GO:GO:0005737',  'GO:GO:0006094', 'GO:GO:0006096',  'GO:GO:0006098'],['Q5FAJ5', 'GO:GO:0004719'],['Q5FAL3', 'GO:GO:0004872',  'GO:GO:0005506', 'GO:GO:0009279',  'GO:GO:0015343'],['Q5FAK9', 'GO:GO:0003700',  'GO:GO:0006351', 'GO:GO:0043565'],['Q5FAK8', 'GO:GO:0004252', 'GO:GO:0070008'],['Q5FAK7', 'GO:GO:0004042',  'GO:GO:0005737', 'GO:GO:0006526'] ,['Q5FAK5', 'GO:GO:0000287', 'GO:GO:0004588', 'GO:GO:0044205'], ['Q5FAI7', 'GO:GO:0006474', 'GO:GO:0008080'],['Q5FAI6', 'GO:GO:0070526'], ['Q5FAI5', 'GO:GO:0009279',  'GO:GO:0016021'], ['Q5FAI4', 'GO:GO:0004332',  'GO:GO:0006096',  'GO:GO:0008270'],['Q5FAI3', 'GO:GO:0003677', 'GO:GO:0005737',  'GO:GO:0006313', 'GO:GO:0007049', 'GO:GO:0007059',  'GO:GO:0009037', 'GO:GO:0051301'], ['Q5FAI2', 'GO:GO:0000287',  'GO:GO:0008661',  'GO:GO:0009228', 'GO:GO:0016114', 'GO:GO:0030976', 'GO:GO:0052865'],['Q5FAI1', 'GO:GO:0005506', 'GO:GO:0005737', 'GO:GO:0006400', 'GO:GO:0016740',  'GO:GO:0051539'], ['Q5FAH9', 'GO:GO:0005737', 'GO:GO:0006782', 'GO:GO:0008483', 'GO:GO:0030170', 'GO:GO:0042286'],['Q5FAH8', 'GO:GO:0003676', 'GO:GO:0005737', 'GO:GO:0016896'],['Q5FAH7', 'GO:GO:0005737', 'GO:GO:0008276'],['Q5FAH6', 'GO:GO:0004075', 'GO:GO:0005524',  'GO:GO:0046872'],['Q5FAH5', 'GO:GO:0003989', 'GO:GO:0006633', 'GO:GO:0009317'], ['Q5FAH3', 'GO:GO:0005737',  'GO:GO:0008616',  'GO:GO:0016740', 'GO:GO:0016853'],['Q5FAH2', 'GO:GO:0004088',  'GO:GO:0005524', 'GO:GO:0006526', 'GO:GO:0044205', 'GO:GO:0046872'],['Q5FAG9', 'GO:GO:0004088',  'GO:GO:0005524',  'GO:GO:0006526', 'GO:GO:0006543', 'GO:GO:0044205', 'GO:GO:0070409'],['Q5FAG5', 'GO:GO:0016491'],['Q5FAG4', 'GO:GO:0003677',  'GO:GO:0003700',  'GO:GO:0005622',  'GO:GO:0006351',  'GO:GO:0010124', 'GO:GO:0045892'],['Q5FAG3', 'GO:GO:0010181',  'GO:GO:0016651', 'GO:GO:0042537',  'GO:GO:0042602',  'GO:GO:0051287'], ['Q5FAG2', 'GO:GO:0009058', 'GO:GO:0016779'],['Q5FAG1', 'GO:GO:0004329', 'GO:GO:0005524', 'GO:GO:0009396', 'GO:GO:0035999'],['Q5FAG0', 'GO:GO:0005524', 'GO:GO:0005525', 'GO:GO:0016887', 'GO:GO:0043022', 'GO:GO:0043023'], ['Q5FAF8', 'GO:GO:0016020',  'GO:GO:0016747'],['Q5FAF7', 'GO:GO:0003723', 'GO:GO:0004831', 'GO:GO:0005524',  'GO:GO:0005737',  'GO:GO:0006437'],['Q5FAF6', 'GO:GO:0003919', 'GO:GO:0008531',  'GO:GO:0009231'], ['Q5FAF5', 'GO:GO:0002161',  'GO:GO:0004822',  'GO:GO:0005524', 'GO:GO:0005737', 'GO:GO:0006428', 'GO:GO:0008270'],['Q5FAF4', 'GO:GO:0009279', 'GO:GO:0015288', 'GO:GO:0016021'], ['Q5FAF3', 'GO:GO:0004190', 'GO:GO:0005886',  'GO:GO:0016021'], ['Q5FAF2', 'GO:GO:0016114', 'GO:GO:0019288', 'GO:GO:0046872',  'GO:GO:0050992', 'GO:GO:0051538', 'GO:GO:0051745'],['Q5FAF1', 'GO:GO:0008967'],[ 'Q5FAE8', 'GO:GO:0003677', 'GO:GO:0003887', 'GO:GO:0005737', 'GO:GO:0006260',  'GO:GO:0008408'],['Q5FAE7', 'GO:GO:0003723'], ['Q5FAE5', 'GO:GO:0003677',  'GO:GO:0003887', 'GO:GO:0006260'], ['Q5FAE4', 'GO:GO:0009042', 'GO:GO:0009058', 'GO:GO:0030170'],['Q5FAE3', 'GO:GO:0009058'],[ 'Q5FAE2', 'GO:GO:0003824', 'GO:GO:0030170'],['Q5FAD9', 'GO:GO:0009058'],['Q5FAD8', 'GO:GO:0000271', 'GO:GO:0016020'], ['Q5FAD7', 'GO:GO:0008270',  'GO:GO:0008703', 'GO:GO:0008835', 'GO:GO:0009231',  'GO:GO:0050661'], ['Q5FAD6', 'GO:GO:0003677', 'GO:GO:0005524', 'GO:GO:0006351',  'GO:GO:0008270', 'GO:GO:0045892'], ['Q5FAD4', 'GO:GO:0003856',  'GO:GO:0005737',  'GO:GO:0009073', 'GO:GO:0009423'],['Q5FAD3', 'GO:GO:0000287',  'GO:GO:0004765', 'GO:GO:0005524',  'GO:GO:0005737', 'GO:GO:0009073',  'GO:GO:0009423',  'GO:GO:0008565'],[ 'Q5FAD2', 'GO:GO:0009279', 'GO:GO:0009297',  'GO:GO:0009306', 'GO:GO:0030420'],['Q5FAC7', 'GO:GO:0005886',  'GO:GO:0008233', 'GO:GO:0008360', 'GO:GO:0008658', 'GO:GO:0009252',  'GO:GO:0016021',  'GO:GO:0016763', 'GO:GO:0046677', 'GO:GO:0071555'], ['Q5FAC6', 'GO:GO:0000287',  'GO:GO:0000917', 'GO:GO:0005525'], ['Q5FAC5', 'GO:GO:0005506', 'GO:GO:0009055',  'GO:GO:0020037', 'GO:GO:0042597', 'Q5FAC3', 'GO:GO:0017004', 'GO:GO:0020037'], ['Q5FAC2', 'GO:GO:0004222', 'GO:GO:0005506', 'GO:GO:0005737', 'GO:GO:0016747', 'GO:GO:0070526'],['Q5FAC1', 'GO:GO:0009244', 'GO:GO:0016021', 'GO:GO:0016746'], ['Q5FAC0', 'GO:GO:0000287', 'GO:GO:0004478' , 'GO:GO:0005524', 'GO:GO:0005737',  'GO:GO:0006556',  'GO:GO:0006730'],[ 'Q5FAB9', 'GO:GO:0009002'],['Q5FAB6', 'GO:GO:0015137', 'GO:GO:0016021'], ['Q5FAB5', 'GO:GO:0000155', 'GO:GO:0005524', 'GO:GO:0016020'],['Q5FAB4', 'GO:GO:0003723', 'GO:GO:0004540', 'GO:GO:0006396'],['Q5FAB3', 'GO:GO:0009055', 'GO:GO:0015035', 'GO:GO:0045454'], ['Q5FAB1', 'GO:GO:0005737', 'GO:GO:0006457', 'GO:GO:0015031', 'GO:GO:0051262'],['Q5FAB0', 'GO:GO:0003676', 'GO:GO:0004003', 'GO:GO:0005524', 'GO:GO:0006281', 'GO:GO:0006310'],['Q5FAA9', 'GO:GO:0003942',  'GO:GO:0005737', 'GO:GO:0006526', 'GO:GO:0051287'],['Q5FAA4', 'GO:GO:0009306', 'GO:GO:0016020', 'GO:GO:0055085'], ['Q5FAA1', 'GO:GO:0005887'],['Q5FA99', 'GO:GO:0005524',  'GO:GO:0008270', 'GO:GO:0008616', 'GO:GO:0016879'],['Q5FA97', 'GO:GO:0008616', 'GO:GO:0016829', 'GO:GO:0046872'], ['Q5FA96', 'GO:GO:0008616',  'GO:GO:0016840', 'GO:GO:0046872', 'GO:GO:0051539'],['Q5FA94', 'GO:GO:0004563', 'GO:GO:0005737', 'GO:GO:0005975', 'GO:GO:0007049', 'GO:GO:0008360', 'GO:GO:0009252',  'GO:GO:0009254', 'GO:GO:0051301', 'GO:GO:0071555'],[ 'Q5FA93', 'GO:GO:0016021'], ['Q5FA91', 'GO:GO:0004252'],['Q5FA90', 'GO:GO:0003677', 'GO:GO:0003906', 'GO:GO:0006284', 'GO:GO:0019104',  'GO:GO:0046872', 'GO:GO:0051539'],[ 'Q5FA87', 'GO:GO:0005354', 'GO:GO:0005355', 'GO:GO:0009276', 'GO:GO:0016021'],['Q5FA86', 'GO:GO:0015297', 'GO:GO:0016021'], ['Q5FA85', 'GO:GO:0016614', 'GO:GO:0030554', 'GO:GO:0050660'],['Q5FA84', 'GO:GO:0004222', 'GO:GO:0004521', 'GO:GO:0005737', 'GO:GO:0006364', 'GO:GO:0008270'],['Q5FA83', 'GO:GO:0004418', 'GO:GO:0006782',  'GO:GO:0018160'], ['Q5FA80', 'GO:GO:0003676', 'GO:GO:0005524', 'GO:GO:0008026'], ['Q5FA78', 'GO:GO:0017150', 'GO:GO:0050660'], ['Q5FA76', 'GO:GO:0000287', 'GO:GO:0003676', 'GO:GO:0006281', 'GO:GO:0006310', 'GO:GO:0008821'],['Q5FA75', 'GO:GO:0009244', 'GO:GO:0016021',  'GO:GO:0016746'], ['Q5FA73', 'GO:GO:0003676'], ['Q5FA72', 'GO:GO:0004177', 'GO:GO:0008237',  'GO:GO:0008270'], ['Q5FA65', 'GO:GO:0009279', 'GO:GO:0016021'],['Q5FA63', 'GO:GO:0007155', 'GO:GO:0030001', 'GO:GO:0046872'], ['Q5FA62', 'GO:GO:0005524', 'GO:GO:0016021', 'GO:GO:0042626'],['Q5FA61', 'GO:GO:0005524', 'GO:GO:0016887'], ['Q5FA60', 'GO:GO:0003735', 'GO:GO:0005840',  'GO:GO:0006412'],['Q5FA59', 'GO:GO:0005737', 'GO:GO:0052906'], ['Q5FA58', 'GO:GO:0005840', 'GO:GO:0006364', 'GO:GO:0042274', 'GO:GO:0043022'],['Q5FA57', 'GO:GO:0003735',  'GO:GO:0005840','GO:GO:0006412'],['Q5FA56', 'GO:GO:0000155',  'GO:GO:0005524', 'GO:GO:0016021'], ['Q5FA55', 'GO:GO:0000160',  'GO:GO:0003677', 'GO:GO:0006351', 'GO:GO:0006355'],['Q5FA52', 'GO:GO:0005737'], ['Q5FA51', 'GO:GO:0005887', 'GO:GO:0008320',  'GO:GO:0033281', 'GO:GO:0043953'],['Q5FA50', 'GO:GO:0005887', 'GO:GO:0008320',  'GO:GO:0009306', 'GO:GO:0033281',  'GO:GO:0043953'],['Q5FA49', 'GO:GO:0005887', 'GO:GO:0008320', 'GO:GO:0009306', 'GO:GO:0033281',  'GO:GO:0043953'],['Q5FA48', 'GO:GO:0003824'],[ 'Q5FA47', 'GO:GO:0000105', 'GO:GO:0004636', 'GO:GO:0005524', 'GO:GO:0005737'],['Q5FA46', 'GO:GO:0008270',  'GO:GO:0016491'], ['Q5FA43', 'GO:GO:0005622',  'GO:GO:0005886', 'GO:GO:0006605',  'GO:GO:0015450', 'GO:GO:0016021',  'GO:GO:0043952', 'GO:GO:0065002'],['Q5FA42', 'GO:GO:0005622',  'GO:GO:0005886', 'GO:GO:0006605', 'GO:GO:0015450', 'GO:GO:0016021', 'GO:GO:0043952',  'GO:GO:0065002'], ['Q5FA41', 'GO:GO:0003735', 'GO:GO:0005840', 'GO:GO:0006412',  'GO:GO:0019843'],['Q5FA40', 'GO:GO:0005524', 'GO:GO:0015417', 'GO:GO:0043190'],[ 'Q5FA39', 'GO:GO:0006810', 'GO:GO:0016021'],[ 'Q5FA38', 'GO:GO:0006810',  'GO:GO:0016021'] ,['Q5FA37', 'GO:GO:0016491'], ['Q5FA36', 'GO:GO:0008519', 'GO:GO:0016020'],['Q5FA35', 'GO:GO:0003723', 'GO:GO:0004386',  'GO:GO:0005524', 'GO:GO:0006353', 'GO:GO:0006355',  'GO:GO:0008186'],[ 'Q5FA34', 'GO:GO:0005524', 'GO:GO:0006090',  'GO:GO:0006094', 'GO:GO:0008986', 'GO:GO:0046872'],['Q5FA32', 'GO:GO:0004674',  'GO:GO:0005524'], ['Q5FA31', 'GO:GO:0016787'], ['Q5FA29', 'GO:GO:0042597', 'GO:GO:0042954'],['Q5FA28', 'GO:GO:0015846', 'GO:GO:0019808', 'GO:GO:0042597'], ['Q5FA27', 'GO:GO:0016740'],['Q5FA25', 'GO:GO:0003924', 'GO:GO:0005525', 'GO:GO:0005737', 'GO:GO:0006449',  'GO:GO:0016149'],['Q5FA24', 'GO:GO:0000105', 'GO:GO:0004635',  'GO:GO:0005737', 'GO:GO:0046872'],['Q5FA23', 'GO:GO:0000105', 'GO:GO:0000107', 'GO:GO:0005737', 'GO:GO:0016829'],['Q5FA22', 'GO:GO:0000105', 'GO:GO:0003949', 'GO:GO:0005737'],['Q5FA21', 'GO:GO:0000105',  'GO:GO:0000107', 'GO:GO:0005737', 'GO:GO:0006541'],['Q5FA20', 'GO:GO:0016407'],['Q5FA19', 'GO:GO:0005524', 'GO:GO:0015408', 'GO:GO:0043190', 'GO:GO:0055072'],['Q5FA18', 'GO:GO:0006810', 'GO:GO:0016021'], ['Q5FA17', 'GO:GO:0005215'], ['Q5FA15', 'GO:GO:0004056', 'GO:GO:0005737', 'GO:GO:0042450'],['Q5FA14', 'GO:GO:0003983', 'GO:GO:0006011',  'GO:GO:0009058'], ['Q5FA13', 'GO:GO:0000166',  'GO:GO:0006163', 'GO:GO:0009143', 'GO:GO:0017111', 'GO:GO:0046872', 'GO:GO:0047429'], ['Q5FA11', 'GO:GO:0000287', 'GO:GO:0004427', 'GO:GO:0005737',  'GO:GO:0006796'],['Q5FA10', 'GO:GO:0006281', 'GO:GO:0008828'], ['Q5FA07', 'GO:GO:0004519'], ['Q5FA04', 'GO:GO:0005886',  'GO:GO:0016021', 'GO:GO:0022820'], ['Q5FA03', 'GO:GO:0008803'],['Q5FA02', 'GO:GO:0019239'],['Q5FA01', 'GO:GO:0009279', 'GO:GO:0015288', 'GO:GO:0016021'], ['Q5FA00', 'GO:GO:0004109', 'GO:GO:0005737', 'GO:GO:0006779', 'GO:GO:0051536'],['Q5F9Z9', 'GO:GO:0003677', 'GO:GO:0003911',  'GO:GO:0006260', 'GO:GO:0006281', 'GO:GO:0046872'], ['Q5F9Z8', 'GO:GO:0000917', 'GO:GO:0005886',  'GO:GO:0016021'],[ 'Q5F9Z7', 'GO:GO:0008745', 'GO:GO:0009253'], ['Q5F9Z5', 'GO:GO:0004798', 'GO:GO:0005524', 'GO:GO:0006233',  'GO:GO:0006235'],['Q5F9Z4', 'GO:GO:0004471', 'GO:GO:0006108',  'GO:GO:0008948', 'GO:GO:0046872',  'GO:GO:0051287'], ['Q5F9Z2', 'GO:GO:0005524',  'GO:GO:0009029', 'GO:GO:0009245']]

#get GO terms from URL
def get_GO_terms(lista):
    nome=[]
    for i in range(len(lista)):
        nome.append([])
        for j in range(1,len(lista[i])):
            GO=str(lista[i][j])
            data = urllib.request.urlopen("http://www.ebi.ac.uk/QuickGO/GTerm?id="+GO[3:13]).read()
            nome[i].append(data.split()[95])
            nome[i].append(data.split()[96])
            nome[i].append(data.split()[97])
            nome[i].append(data.split()[98])
            nome[i].append(data.split()[99])
            nome[i].append(data.split()[100])
    return nome

G=[[b'<title>GO:0003688', b'DNA', b'replication', b'origin', b'binding</title>', b'<script', b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006270', b'DNA', b'replication', b'initiation</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0006275', b'regulation', b'of', b'DNA', b'replication</title>', b'<script'], [b'<title>GO:0003677', b'DNA', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0003887', b'DNA-directed', b'DNA', b'polymerase', b'activity</title>', b'<script', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006260', b'DNA', b'replication</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0008408', b"3'-5'", b'exonuclease', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0009360', b'DNA', b'polymerase', b'III', b'complex</title>', b'<script'], [b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0006799', b'polyphosphate', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0008976', b'polyphosphate', b'kinase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0009358', b'polyphosphate', b'kinase', b'complex</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0002161', b'aminoacyl-tRNA', b'editing', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0004823', b'leucine-tRNA', b'ligase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006429', b'leucyl-tRNA', b'aminoacylation</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0043565', b'sequence-specific', b'DNA', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0009306', b'protein', b'secretion</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0015450', b'P-P-bond-hydrolysis-driven', b'protein', b'transmembrane', b'transporter', b'activity</title>', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script'], [b'<title>GO:0004807', b'triose-phosphate', b'isomerase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006094', b'gluconeogenesis</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006096', b'glycolytic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0006098', b'pentose-phosphate', b'shunt</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0004719', b'protein-L-isoaspartate', b'(D-aspartate)', b'O-methyltransferase', b'activity</title>', b'<script'], [b'<title>GO:0004872', b'receptor', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0005506', b'iron', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0009279', b'cell', b'outer', b'membrane</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0015343', b'siderophore', b'transmembrane', b'transporter', b'activity</title>', b'<script'], [b'<title>GO:0003700', b'sequence-specific', b'DNA', b'binding', b'transcription', b'factor', b'<title>GO:0006351', b'transcription,', b'DNA-templated</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0043565', b'sequence-specific', b'DNA', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0004252', b'serine-type', b'endopeptidase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0070008', b'serine-type', b'exopeptidase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0004042', b'acetyl-CoA:L-glutamate', b'N-acetyltransferase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006526', b'arginine', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0000287', b'magnesium', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0004588', b'orotate', b'phosphoribosyltransferase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0044205', b"'de", b"novo'", b'UMP', b'biosynthetic', b'process</title>'], [b'<title>GO:0006474', b'N-terminal', b'protein', b'amino', b'acid', b'acetylation</title>', b'<title>GO:0008080', b'N-acetyltransferase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0070526', b'threonylcarbamoyladenosine', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0009279', b'cell', b'outer', b'membrane</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script'], [b'<title>GO:0004332', b'fructose-bisphosphate', b'aldolase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0006096', b'glycolytic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0008270', b'zinc', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0003677', b'DNA', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006313', b'transposition,', b'DNA-mediated</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0007049', b'cell', b'cycle</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0007059', b'chromosome', b'segregation</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0009037', b'tyrosine-based', b'site-specific', b'recombinase', b'activity</title>', b'<script', b'<title>GO:0051301', b'cell', b'division</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0000287', b'magnesium', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0008661', b'1-deoxy-D-xylulose-5-phosphate', b'synthase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0009228', b'thiamine', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0016114', b'terpenoid', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0030976', b'thiamine', b'pyrophosphate', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0052865', b'1-deoxy-D-xylulose', b'5-phosphate', b'biosynthetic', b'process</title>', b'<script'], [b'<title>GO:0005506', b'iron', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006400', b'tRNA', b'modification</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0016740', b'transferase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0051539', b'4', b'iron,', b'4', b'sulfur', b'cluster'], [b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006782', b'protoporphyrinogen', b'IX', b'biosynthetic', b'process</title>', b'<script', b'<title>GO:0008483', b'transaminase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0030170', b'pyridoxal', b'phosphate', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0042286', b'glutamate-1-semialdehyde', b'2,1-aminomutase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0003676', b'nucleic', b'acid', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0016896', b'exoribonuclease', b'activity,', b'producing', b"5'-phosphomonoesters</title>", b'<script'], [b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0008276', b'protein', b'methyltransferase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0004075', b'biotin', b'carboxylase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0046872', b'metal', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0003989', b'acetyl-CoA', b'carboxylase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0006633', b'fatty', b'acid', b'biosynthetic', b'process</title>', b'<script', b'<title>GO:0009317', b'acetyl-CoA', b'carboxylase', b'complex</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0008616', b'queuosine', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0016740', b'transferase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0016853', b'isomerase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0004088', b'carbamoyl-phosphate', b'synthase', b'(glutamine-hydrolyzing)', b'activity</title>', b'<script', b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0006526', b'arginine', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0044205', b"'de", b"novo'", b'UMP', b'biosynthetic', b'process</title>', b'<title>GO:0046872', b'metal', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0004088', b'carbamoyl-phosphate', b'synthase', b'(glutamine-hydrolyzing)', b'activity</title>', b'<script', b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0006526', b'arginine', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0006543', b'glutamine', b'catabolic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0044205', b"'de", b"novo'", b'UMP', b'biosynthetic', b'process</title>', b'<title>GO:0070409', b'carbamoyl', b'phosphate', b'biosynthetic', b'process</title>', b'<script'], [b'<title>GO:0016491', b'oxidoreductase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0003677', b'DNA', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0003700', b'sequence-specific', b'DNA', b'binding', b'transcription', b'factor', b'<title>GO:0005622', b'intracellular</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006351', b'transcription,', b'DNA-templated</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0010124', b'phenylacetate', b'catabolic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0045892', b'negative', b'regulation', b'of', b'transcription,', b'DNA-templated</title>'], [b'<title>GO:0010181', b'FMN', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0016651', b'oxidoreductase', b'activity,', b'acting', b'on', b'NAD(P)H</title>', b'<title>GO:0042537', b'benzene-containing', b'compound', b'metabolic', b'process</title>', b'<script', b'<title>GO:0042602', b'riboflavin', b'reductase', b'(NADPH)', b'activity</title>', b'<script', b'<title>GO:0051287', b'NAD', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0009058', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0016779', b'nucleotidyltransferase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0004329', b'formate-tetrahydrofolate', b'ligase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0009396', b'folic', b'acid-containing', b'compound', b'biosynthetic', b'process</title>', b'<title>GO:0035999', b'tetrahydrofolate', b'interconversion</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0005525', b'GTP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0016887', b'ATPase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0043022', b'ribosome', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0043023', b'ribosomal', b'large', b'subunit', b'binding</title>', b'<script'], [b'<title>GO:0016020', b'membrane</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0016747', b'transferase', b'activity,', b'transferring', b'acyl', b'groups'], [b'<title>GO:0003723', b'RNA', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0004831', b'tyrosine-tRNA', b'ligase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006437', b'tyrosyl-tRNA', b'aminoacylation</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0003919', b'FMN', b'adenylyltransferase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0008531', b'riboflavin', b'kinase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0009231', b'riboflavin', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0002161', b'aminoacyl-tRNA', b'editing', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0004822', b'isoleucine-tRNA', b'ligase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006428', b'isoleucyl-tRNA', b'aminoacylation</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0008270', b'zinc', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0009279', b'cell', b'outer', b'membrane</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0015288', b'porin', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script'], [b'<title>GO:0004190', b'aspartic-type', b'endopeptidase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005886', b'plasma', b'membrane</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script'], [b'<title>GO:0016114', b'terpenoid', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0019288', b'isopentenyl', b'diphosphate', b'biosynthetic', b'process,', b'methylerythritol', b'<title>GO:0046872', b'metal', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0050992', b'dimethylallyl', b'diphosphate', b'biosynthetic', b'process</title>', b'<script', b'<title>GO:0051538', b'3', b'iron,', b'4', b'sulfur', b'cluster', b'<title>GO:0051745', b'4-hydroxy-3-methylbut-2-en-1-yl', b'diphosphate', b'reductase', b'activity</title>', b'<script'], [b'<title>GO:0008967', b'phosphoglycolate', b'phosphatase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0003677', b'DNA', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0003887', b'DNA-directed', b'DNA', b'polymerase', b'activity</title>', b'<script', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006260', b'DNA', b'replication</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0008408', b"3'-5'", b'exonuclease', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0003723', b'RNA', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0003677', b'DNA', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0003887', b'DNA-directed', b'DNA', b'polymerase', b'activity</title>', b'<script', b'<title>GO:0006260', b'DNA', b'replication</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0009042', b'valine-pyruvate', b'transaminase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0009058', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0030170', b'pyridoxal', b'phosphate', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0009058', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0003824', b'catalytic', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0030170', b'pyridoxal', b'phosphate', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0009058', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0000271', b'polysaccharide', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0016020', b'membrane</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function'], [b'<title>GO:0008270', b'zinc', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0008703', b'5-amino-6-(5-phosphoribosylamino)uracil', b'reductase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0008835', b'diaminohydroxyphosphoribosylaminopyrimidine', b'deaminase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0009231', b'riboflavin', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0050661', b'NADP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0003677', b'DNA', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0006351', b'transcription,', b'DNA-templated</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0008270', b'zinc', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0045892', b'negative', b'regulation', b'of', b'transcription,', b'DNA-templated</title>'], [b'<title>GO:0003856', b'3-dehydroquinate', b'synthase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0009073', b'aromatic', b'amino', b'acid', b'family', b'biosynthetic', b'<title>GO:0009423', b'chorismate', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0000287', b'magnesium', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0004765', b'shikimate', b'kinase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0009073', b'aromatic', b'amino', b'acid', b'family', b'biosynthetic', b'<title>GO:0009423', b'chorismate', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0008565', b'protein', b'transporter', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0009279', b'cell', b'outer', b'membrane</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0009297', b'pilus', b'assembly</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0009306', b'protein', b'secretion</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0030420', b'establishment', b'of', b'competence', b'for', b'transformation</title>'], [b'<title>GO:0005886', b'plasma', b'membrane</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0008233', b'peptidase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0008360', b'regulation', b'of', b'cell', b'shape</title>', b'<script', b'<title>GO:0008658', b'penicillin', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0009252', b'peptidoglycan', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script', b'<title>GO:0016763', b'transferase', b'activity,', b'transferring', b'pentosyl', b'groups</title>', b'<title>GO:0046677', b'response', b'to', b'antibiotic</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0071555', b'cell', b'wall', b'organization</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0000287', b'magnesium', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0000917', b'barrier', b'septum', b'assembly</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005525', b'GTP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0005506', b'iron', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0009055', b'electron', b'carrier', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0020037', b'heme', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0042597', b'periplasmic', b'space</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO', b'Term', b'not', b'found</title>', b'<script>', b'function', b'<title>GO:0017004', b'cytochrome', b'complex', b'assembly</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0020037', b'heme', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0004222', b'metalloendopeptidase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0005506', b'iron', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0016747', b'transferase', b'activity,', b'transferring', b'acyl', b'groups', b'<title>GO:0070526', b'threonylcarbamoyladenosine', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0009244', b'lipopolysaccharide', b'core', b'region', b'biosynthetic', b'process</title>', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script', b'<title>GO:0016746', b'transferase', b'activity,', b'transferring', b'acyl', b'groups</title>'], [b'<title>GO:0000287', b'magnesium', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0004478', b'methionine', b'adenosyltransferase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006556', b'S-adenosylmethionine', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0006730', b'one-carbon', b'metabolic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0009002', b'serine-type', b'D-Ala-D-Ala', b'carboxypeptidase', b'activity</title>', b'<script'], [b'<title>GO:0015137', b'citrate', b'transmembrane', b'transporter', b'activity</title>', b'<script', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script'], [b'<title>GO:0000155', b'phosphorelay', b'sensor', b'kinase', b'activity</title>', b'<script', b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0016020', b'membrane</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function'], [b'<title>GO:0003723', b'RNA', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0004540', b'ribonuclease', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0006396', b'RNA', b'processing</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0009055', b'electron', b'carrier', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0015035', b'protein', b'disulfide', b'oxidoreductase', b'activity</title>', b'<script', b'<title>GO:0045454', b'cell', b'redox', b'homeostasis</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006457', b'protein', b'folding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0015031', b'protein', b'transport</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0051262', b'protein', b'tetramerization</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0003676', b'nucleic', b'acid', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0004003', b'ATP-dependent', b'DNA', b'helicase', b'activity</title>', b'<script', b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0006281', b'DNA', b'repair</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0006310', b'DNA', b'recombination</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0003942', b'N-acetyl-gamma-glutamyl-phosphate', b'reductase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006526', b'arginine', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0051287', b'NAD', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0009306', b'protein', b'secretion</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0016020', b'membrane</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0055085', b'transmembrane', b'transport</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0005887', b'integral', b'component', b'of', b'plasma', b'membrane</title>'], [b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0008270', b'zinc', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0008616', b'queuosine', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0016879', b'ligase', b'activity,', b'forming', b'carbon-nitrogen', b'bonds</title>'], [b'<title>GO:0008616', b'queuosine', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0016829', b'lyase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0046872', b'metal', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0008616', b'queuosine', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0016840', b'carbon-nitrogen', b'lyase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0046872', b'metal', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0051539', b'4', b'iron,', b'4', b'sulfur', b'cluster'], [b'<title>GO:0004563', b'beta-N-acetylhexosaminidase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0005975', b'carbohydrate', b'metabolic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0007049', b'cell', b'cycle</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0008360', b'regulation', b'of', b'cell', b'shape</title>', b'<script', b'<title>GO:0009252', b'peptidoglycan', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0009254', b'peptidoglycan', b'turnover</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0051301', b'cell', b'division</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0071555', b'cell', b'wall', b'organization</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script'], [b'<title>GO:0004252', b'serine-type', b'endopeptidase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0003677', b'DNA', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0003906', b'DNA-(apurinic', b'or', b'apyrimidinic', b'site)', b'lyase', b'<title>GO:0006284', b'base-excision', b'repair</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0019104', b'DNA', b'N-glycosylase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0046872', b'metal', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0051539', b'4', b'iron,', b'4', b'sulfur', b'cluster'], [b'<title>GO:0005354', b'galactose', b'transmembrane', b'transporter', b'activity</title>', b'<script', b'<title>GO:0005355', b'glucose', b'transmembrane', b'transporter', b'activity</title>', b'<script', b'<title>GO:0009276', b'Gram-negative-bacterium-type', b'cell', b'wall</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script'], [b'<title>GO:0015297', b'antiporter', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script'], [b'<title>GO:0016614', b'oxidoreductase', b'activity,', b'acting', b'on', b'CH-OH', b'<title>GO:0030554', b'adenyl', b'nucleotide', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0050660', b'flavin', b'adenine', b'dinucleotide', b'binding</title>', b'<script'], [b'<title>GO:0004222', b'metalloendopeptidase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0004521', b'endoribonuclease', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006364', b'rRNA', b'processing</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0008270', b'zinc', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0004418', b'hydroxymethylbilane', b'synthase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0006782', b'protoporphyrinogen', b'IX', b'biosynthetic', b'process</title>', b'<script', b'<title>GO:0018160', b'peptidyl-pyrromethane', b'cofactor', b'linkage</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0003676', b'nucleic', b'acid', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0008026', b'ATP-dependent', b'helicase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0017150', b'tRNA', b'dihydrouridine', b'synthase', b'activity</title>', b'<script', b'<title>GO:0050660', b'flavin', b'adenine', b'dinucleotide', b'binding</title>', b'<script'], [b'<title>GO:0000287', b'magnesium', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0003676', b'nucleic', b'acid', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0006281', b'DNA', b'repair</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0006310', b'DNA', b'recombination</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0008821', b'crossover', b'junction', b'endodeoxyribonuclease', b'activity</title>', b'<script'], [b'<title>GO:0009244', b'lipopolysaccharide', b'core', b'region', b'biosynthetic', b'process</title>', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script', b'<title>GO:0016746', b'transferase', b'activity,', b'transferring', b'acyl', b'groups</title>'], [b'<title>GO:0003676', b'nucleic', b'acid', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0004177', b'aminopeptidase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0008237', b'metallopeptidase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0008270', b'zinc', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0009279', b'cell', b'outer', b'membrane</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script'], [b'<title>GO:0007155', b'cell', b'adhesion</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0030001', b'metal', b'ion', b'transport</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0046872', b'metal', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script', b'<title>GO:0042626', b'ATPase', b'activity,', b'coupled', b'to', b'transmembrane'], [b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0016887', b'ATPase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0003735', b'structural', b'constituent', b'of', b'ribosome</title>', b'<script', b'<title>GO:0005840', b'ribosome</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006412', b'translation</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function'], [b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0052906', b'tRNA', b'(guanine(37)-N(1))-methyltransferase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0005840', b'ribosome</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006364', b'rRNA', b'processing</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0042274', b'ribosomal', b'small', b'subunit', b'biogenesis</title>', b'<script', b'<title>GO:0043022', b'ribosome', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0003735', b'structural', b'constituent', b'of', b'ribosome</title>', b'<script', b'<title>GO:0005840', b'ribosome</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006412', b'translation</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function'], [b'<title>GO:0000155', b'phosphorelay', b'sensor', b'kinase', b'activity</title>', b'<script', b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script'], [b'<title>GO:0000160', b'phosphorelay', b'signal', b'transduction', b'system</title>', b'<script', b'<title>GO:0003677', b'DNA', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0006351', b'transcription,', b'DNA-templated</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0006355', b'regulation', b'of', b'transcription,', b'DNA-templated</title>', b'<script'], [b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function'], [b'<title>GO:0005887', b'integral', b'component', b'of', b'plasma', b'membrane</title>', b'<title>GO:0008320', b'protein', b'transmembrane', b'transporter', b'activity</title>', b'<script', b'<title>GO:0033281', b'TAT', b'protein', b'transport', b'complex</title>', b'<script', b'<title>GO:0043953', b'protein', b'transport', b'by', b'the', b'Tat'], [b'<title>GO:0005887', b'integral', b'component', b'of', b'plasma', b'membrane</title>', b'<title>GO:0008320', b'protein', b'transmembrane', b'transporter', b'activity</title>', b'<script', b'<title>GO:0009306', b'protein', b'secretion</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0033281', b'TAT', b'protein', b'transport', b'complex</title>', b'<script', b'<title>GO:0043953', b'protein', b'transport', b'by', b'the', b'Tat'], [b'<title>GO:0005887', b'integral', b'component', b'of', b'plasma', b'membrane</title>', b'<title>GO:0008320', b'protein', b'transmembrane', b'transporter', b'activity</title>', b'<script', b'<title>GO:0009306', b'protein', b'secretion</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0033281', b'TAT', b'protein', b'transport', b'complex</title>', b'<script', b'<title>GO:0043953', b'protein', b'transport', b'by', b'the', b'Tat'], [b'<title>GO:0003824', b'catalytic', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0000105', b'histidine', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0004636', b'phosphoribosyl-ATP', b'diphosphatase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function'], [b'<title>GO:0008270', b'zinc', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0016491', b'oxidoreductase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0005622', b'intracellular</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0005886', b'plasma', b'membrane</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0006605', b'protein', b'targeting</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0015450', b'P-P-bond-hydrolysis-driven', b'protein', b'transmembrane', b'transporter', b'activity</title>', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script', b'<title>GO:0043952', b'protein', b'transport', b'by', b'the', b'Sec', b'<title>GO:0065002', b'intracellular', b'protein', b'transmembrane', b'transport</title>', b'<script'], [b'<title>GO:0005622', b'intracellular</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0005886', b'plasma', b'membrane</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0006605', b'protein', b'targeting</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0015450', b'P-P-bond-hydrolysis-driven', b'protein', b'transmembrane', b'transporter', b'activity</title>', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script', b'<title>GO:0043952', b'protein', b'transport', b'by', b'the', b'Sec', b'<title>GO:0065002', b'intracellular', b'protein', b'transmembrane', b'transport</title>', b'<script'], [b'<title>GO:0003735', b'structural', b'constituent', b'of', b'ribosome</title>', b'<script', b'<title>GO:0005840', b'ribosome</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006412', b'translation</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0019843', b'rRNA', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0015417', b'polyamine-transporting', b'ATPase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0043190', b'ATP-binding', b'cassette', b'(ABC)', b'transporter', b'complex</title>'], [b'<title>GO:0006810', b'transport</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script'], [b'<title>GO:0006810', b'transport</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script'], [b'<title>GO:0016491', b'oxidoreductase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0008519', b'ammonium', b'transmembrane', b'transporter', b'activity</title>', b'<script', b'<title>GO:0016020', b'membrane</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function'], [b'<title>GO:0003723', b'RNA', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0004386', b'helicase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0006353', b'DNA-templated', b'transcription,', b'termination</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0006355', b'regulation', b'of', b'transcription,', b'DNA-templated</title>', b'<script', b'<title>GO:0008186', b'RNA-dependent', b'ATPase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0006090', b'pyruvate', b'metabolic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0006094', b'gluconeogenesis</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0008986', b'pyruvate,', b'water', b'dikinase', b'activity</title>', b'<script', b'<title>GO:0046872', b'metal', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0004674', b'protein', b'serine/threonine', b'kinase', b'activity</title>', b'<script', b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0016787', b'hydrolase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0042597', b'periplasmic', b'space</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0042954', b'lipoprotein', b'transporter', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0015846', b'polyamine', b'transport</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0019808', b'polyamine', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0042597', b'periplasmic', b'space</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0016740', b'transferase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0003924', b'GTPase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0005525', b'GTP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006449', b'regulation', b'of', b'translational', b'termination</title>', b'<script', b'<title>GO:0016149', b'translation', b'release', b'factor', b'activity,', b'codon'], [b'<title>GO:0000105', b'histidine', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0004635', b'phosphoribosyl-AMP', b'cyclohydrolase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0046872', b'metal', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0000105', b'histidine', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0000107', b'imidazoleglycerol-phosphate', b'synthase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0016829', b'lyase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0000105', b'histidine', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0003949', b'1-(5-phosphoribosyl)-5-[(5-phosphoribosylamino)methylideneamino]imidazole-4-carboxamide', b'isomerase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function'], [b'<title>GO:0000105', b'histidine', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0000107', b'imidazoleglycerol-phosphate', b'synthase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006541', b'glutamine', b'metabolic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0016407', b'acetyltransferase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0015408', b'ferric-transporting', b'ATPase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0043190', b'ATP-binding', b'cassette', b'(ABC)', b'transporter', b'complex</title>', b'<title>GO:0055072', b'iron', b'ion', b'homeostasis</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0006810', b'transport</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script'], [b'<title>GO:0005215', b'transporter', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0004056', b'argininosuccinate', b'lyase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0042450', b'arginine', b'biosynthetic', b'process', b'via', b'ornithine</title>'], [b'<title>GO:0003983', b'UTP:glucose-1-phosphate', b'uridylyltransferase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0006011', b'UDP-glucose', b'metabolic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0009058', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0000166', b'nucleotide', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0006163', b'purine', b'nucleotide', b'metabolic', b'process</title>', b'<script', b'<title>GO:0009143', b'nucleoside', b'triphosphate', b'catabolic', b'process</title>', b'<script', b'<title>GO:0017111', b'nucleoside-triphosphatase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0046872', b'metal', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0047429', b'nucleoside-triphosphate', b'diphosphatase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0000287', b'magnesium', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0004427', b'inorganic', b'diphosphatase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006796', b'phosphate-containing', b'compound', b'metabolic', b'process</title>', b'<script'], [b'<title>GO:0006281', b'DNA', b'repair</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0008828', b'dATP', b'pyrophosphohydrolase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0004519', b'endonuclease', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0005886', b'plasma', b'membrane</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script', b'<title>GO:0022820', b'potassium', b'ion', b'symporter', b'activity</title>', b'<script'], [b'<title>GO:0008803', b"bis(5'-nucleosyl)-tetraphosphatase", b'(symmetrical)', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0019239', b'deaminase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0009279', b'cell', b'outer', b'membrane</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0015288', b'porin', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script'], [b'<title>GO:0004109', b'coproporphyrinogen', b'oxidase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005737', b'cytoplasm</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'function', b'<title>GO:0006779', b'porphyrin-containing', b'compound', b'biosynthetic', b'process</title>', b'<script', b'<title>GO:0051536', b'iron-sulfur', b'cluster', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0003677', b'DNA', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0003911', b'DNA', b'ligase', b'(NAD+)', b'activity</title>', b'<script', b'<title>GO:0006260', b'DNA', b'replication</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0006281', b'DNA', b'repair</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0046872', b'metal', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0000917', b'barrier', b'septum', b'assembly</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005886', b'plasma', b'membrane</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0016021', b'integral', b'component', b'of', b'membrane</title>', b'<script'], [b'<title>GO:0008745', b'N-acetylmuramoyl-L-alanine', b'amidase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0009253', b'peptidoglycan', b'catabolic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0004798', b'thymidylate', b'kinase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0006233', b'dTDP', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0006235', b'dTTP', b'biosynthetic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>'], [b'<title>GO:0004471', b'malate', b'dehydrogenase', b'(decarboxylating)', b'(NAD+)', b'activity</title>', b'<title>GO:0006108', b'malate', b'metabolic', b'process</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0008948', b'oxaloacetate', b'decarboxylase', b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0046872', b'metal', b'ion', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0051287', b'NAD', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>'], [b'<title>GO:0005524', b'ATP', b'binding</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<script>', b'<title>GO:0009029', b'tetraacyldisaccharide', b"4'-kinase", b'activity</title>', b'<script', b'src="./js/quickgoAnnotation.js"></script>', b'<title>GO:0009245', b'lipid', b'A', b'biosynthetic', b'process</title>', b'<script']]


#viewing GO and genes
def which_GO(lista,G,pos):
   return lista[pos][0],G[pos]
                  
def blastanaliserwithnote():
    blast=[]    
    for file in os.listdir("../res/blast_with_note"):
        if file.endswith(".xml"):
            blast.append(file)
    E_VALUE_THRESH = 0.05
    lista=[]
    for i in range(len(blast)):
        lista.append([])
        lista[i].append(blast[i])
        result_handle = open("../res/blast_with_note/"+blast[i])
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
    dst = "../res/blast_with_note/nomatch/"#destination folder
    shutil.move(src, dst)
                     
    save_file = open('matches.txt', "w")
    for i in range(len(lista)):              
        if len(lista[i])>2:
            save_file.write(str(lista[i])+'\n')
    save_file.close()
    #        #moving the file to another directory
    path=os.getcwd()
    src = path+"/"+'matches.txt' #source folder
    dst = "../res/blast_with_note/match/"#destination folder
    shutil.move(src, dst)

def uniprotallhitswithnote():
    handle = open("../res/blast_with_note/match/allhits/allhits.txt").readlines()
    
    
    
    for n in range(len (handle)):       
            x=handle[n].split()
            for k in range(len(x)-1,len(x)):
                e=x[k]
                limpo=e.replace("']","")
                gi=(limpo[1:])
                if not os.path.exists("../res/blast_with_note/match/function/teste/"+gi):
                    os.makedirs("../res/blast_with_note/match/function/teste/"+gi)
            for j in range(len(x)-1):
                m=x[j]
                q=m.replace("[" ,"")
                protein=q[1:7]
                site = urllib.request.urlopen("http://www.uniprot.org/uniprot/"+protein+".txt")
                data = site.readlines()
                file = open("../res/blast_with_note/match/function/teste/"+protein+'.txt',"wb") #open file in binary mode
                file.writelines(data)
                file.close()
                try:
                    src = "../res/blast_with_note/match/function/teste/"+protein+'.txt' #source folder
                    dst = "../res/blast_with_note/match/function/teste/"+gi #destination folder
                    shutil.move(src, dst)
                except:
                        pass
def allfunctionswithnote():
    handle = open("../res/blast_with_note/match/allhits/allhits.txt").readlines()
    for n in range(len (handle)):       
            x=handle[n].split()
            for k in range(len(x)-1,len(x)):
                e=x[k]
                limpo=e.replace("']","")
                ginote=(limpo[1:])
                blast=[]    
                for file in os.listdir("../res/blast_with_note/match/function/teste/"+ginote):
                    if file.endswith(".txt"):
                           blast.append(file)
                    lista=[]
                    for j in range(len(blast)):
                        nome=blast[j]
                        gi = nome.replace(".txt","")
                        first = '-!- FUNCTION:'
                        last = 'CC'
                        file = open("../res/blast_with_note/match/function/teste/"+ginote+'/'+gi+'.txt').read()
                        data = file.replace("\n", " ") 
                        try:
                            start = data.rindex( first ) + len( first )
                            end = data.rindex(last, start)
                            novo= data[start:end] 
                            lista.append('Gi: '+gi+'  '+'Possivel função:  '+novo)
                        except:
                            pass   
                        
                    file = open("../res/blast_with_note/match/allhits/funcao_all_hits/"+ginote+".txt",'w')       
                    for i in range(len(lista)):
                        file.write("%s\n" % lista[i])
                    file.close()

    
def menu(record):
    ans=True
    while ans:
        print("""
    1.Accessing info from one gene
    2.without note (returns proteinID)
    3.translation (needs proteinID)
    4.tRNA
    5.pseudogenes
    6.hypothetical proteins
    7.Article with a gene reference
    8.Running protein Blast (needs GI number)
    9.Parsing blast
    10.Uniprot_ID 
    11.Uniprot Identifier, Definition, Subcellular location
    12.Get gi from protein without note
    13.tabela
    14.Get note from protein without note
    15.blastanaliser
    16.pega no primeiro hit do blast
    17.GO numbers
    18.function to see possible function of proteins that didnt have note:
    19.return all hits from blast
    20.Go to uniprot and download information for all hits
    21.Get all information for all hits
    22.GI numbers from genes with note
    23.Get all GO terms
    24.Viewing each position/gene GO terms
    25.Blast gi with note  
    26.Blast analiser with note
    28.return all hits from blast with notes
    29.Go to uniprot and download information for all hits with notes
    30.Get all information for all hits with notes
    31.CDD
    60.Exit
    """)
        ans=input("Choose an option? ")
        if ans=="1":
            locus=locus_tag(record)            
            fa=info(record,locus)
            nr=int(input("Qual gene?"))
            print(aceder(fa,nr))
#            print(locus_tag(record))
#            print(genes_names(record))
        elif ans=="2":
            print(len(without_note(record)))
        elif ans=="3":
            prot_ID=str(input("Protein ID: "))
            print(translation(record,prot_ID))
        elif ans=="4":
            loc=tRNA(record)
            print(info_tRNA(record,loc))
        elif ans=="5":
            locus=pseudogenes(record)
            print(info_pseudogenes(record,locus))
        elif ans=="6":
            print(hypoth_proteins(record))
        elif ans=="7":
            gene=str(input("Gene name: "))
            print(DB_pubmed(gene))
        elif ans=="8":
            file=str(input("Qual o nome a colocar no ficheiro? "))+".xml"
            while os.path.isfile("../res/"+file):
                 file=str(input("Qual o nome a colocar no ficheiro? "))+".xml"
            else:False
            GI=str(input("Qual o GI number da sequencia? "))
            blast(GI,file)
        elif ans=="9":
            file=str(input("Qual o nome do ficheiro? "))+".xml"
            parse_blast(file)
        elif ans=="10":
            #print(uniprot_ID(record))
            print(ID)
        elif ans=="11":
            print(info_uniprot())
            print(more_info_uniprot())
        elif ans=="12":
            print(giwithout_note(record))
        elif ans=="13":
            locus=locus_tag(record)            
            fa=info(record,locus)
            print(tabela(fa,locus))
            locus_tRNA=tRNA(record)
            lista_tRNA=info_tRNA(record,locus_tRNA)
            print(tabela_tRNA(lista_tRNA,locus_tRNA))
            locus_pseudo=pseudogenes(record)
            lista_pseudo=info_pseudogenes(record,locus_pseudo)
            print(tabela_pseudogenes(lista_pseudo,locus_pseudo))
            print(tabela_uniprot())
            print(tabela_uniprot2())
        elif ans=="14":
            blastnote(filename)
        elif ans=="15":
            blastanaliser()
        elif ans=="16":
            gimatch()
        elif ans=="17":
            print(GO())
            
        elif ans=="18":
           file = open("../res/blast_without_note/match/allhits/funcao.txt",'w')
           lista=getfunction()        
           for i in range(len(lista)):
               file.write("%s\n" % lista[i])
#           print(getfunction())

        elif ans=="19":
           file = open("../res/blast_without_note/match/allhits/allhits.txt",'w')
           lista=allhits()         
           for i in range(len(lista)):
               file.write("%s\n" % lista[i])
           file.close()
        elif ans=="20":
            uniprotallhits()
            
        elif ans=="21":    
            allfunctions()
        elif ans=="22":    
            locus=locus_tag(record)
            print(GInumbers(record,locus))

        elif ans=="25":   
            blastwithnote()

        elif ans=="23":    
            print(get_GO_terms(lista_GO))
        elif ans=="24":
            pos=int(input("position?from 0 to 146: "))
            print(which_GO(lista_GO,G,pos))
        elif ans=="26":
            blastanaliserwithnote()
        
        elif ans=="28":
           file = open("../res/blast_with_note/match/allhits/allhits.txt",'w')
           lista=allhitswithnote()         
           for i in range(len(lista)):
               file.write("%s\n" % lista[i])
           file.close()
        elif ans=="29":
            uniprotallhitswithnote()
            
        elif ans=="30":    
            allfunctionswithnote()
        
        elif ans=="31":
            #locus=locus_tag(record)
            #GI=GI_number(record,locus)
            #print(CDD(GI))
            GI=[['GI:59800474', 'CDD:256535', 'CDD:273037', 'CDD:99707', 'CDD:99707', 'CDD:99707', 'CDD:99707', 'CDD:99707', 'CDD:119330', 'CDD:119330'], ['GI:59800475', 'CDD:235541', 'CDD:238082', 'CDD:238082', 'CDD:238082', 'CDD:238082', 'CDD:238082'], ['GI:59800476', 'CDD:235469', 'CDD:257492', 'CDD:251338', 'CDD:197262', 'CDD:197262', 'CDD:197262', 'CDD:197262', 'CDD:197265', 'CDD:197265', 'CDD:197265', 'CDD:197265'], ['GI:59800477'], ['GI:59800478', 'CDD:234743', 'CDD:173906', 'CDD:173906', 'CDD:257918', 'CDD:275460', 'CDD:173912', 'CDD:173912', 'CDD:153412', 'CDD:153412'], ['GI:745998703'], ['GI:59800483', 'CDD:225495', 'CDD:238045', 'CDD:238045', 'CDD:238045', 'CDD:238045'], ['GI:59800485', 'CDD:252199'], ['GI:59800486', 'CDD:238190', 'CDD:238190', 'CDD:238190', 'CDD:238190'], ['GI:59800487', 'CDD:224867'], ['GI:59800488', 'CDD:225316', 'CDD:100107', 'CDD:100107'], ['GI:59800489', 'CDD:238786', 'CDD:238786'], ['GI:59800490', 'CDD:227114', 'CDD:238657', 'CDD:238657', 'CDD:238657'], ['GI:59800491', 'CDD:223734'], ['GI:59800492', 'CDD:223687', 'CDD:238560', 'CDD:238560'], ['GI:59800494', 'CDD:225117', 'CDD:249644'], ['GI:59800495', 'CDD:224422', 'CDD:249771'], ['GI:59800496', 'CDD:273856', 'CDD:239770', 'CDD:239770', 'CDD:239770', 'CDD:239770', 'CDD:173926', 'CDD:173926'], ['GI:59800497', 'CDD:264474'], ['GI:59800498', 'CDD:206754', 'CDD:206754'], ['GI:59800499', 'CDD:275600', 'CDD:198424'], ['GI:59800500', 'CDD:223532', 'CDD:173926', 'CDD:173926', 'CDD:276195'], ['GI:59800501', 'CDD:224135', 'CDD:250152'], ['GI:59800502', 'CDD:255039'], ['GI:59800503', 'CDD:238477', 'CDD:238477', 'CDD:238477', 'CDD:238477', 'CDD:238477'], ['GI:59800504', 'CDD:227307', 'CDD:145844', 'CDD:271179', 'CDD:271179'], ['GI:59800505', 'CDD:235470', 'CDD:238965', 'CDD:238965', 'CDD:132916', 'CDD:132916', 'CDD:132916', 'CDD:132916', 'CDD:251527'], ['GI:59800506', 'CDD:237673', 'CDD:250224', 'CDD:100105', 'CDD:100105', 'CDD:110897'], ['GI:59800508', 'CDD:234607', 'CDD:99735', 'CDD:99735', 'CDD:99735', 'CDD:99735'], ['GI:59800509', 'CDD:99838', 'CDD:99838', 'CDD:99838', 'CDD:99838', 'CDD:99838'], ['GI:59800510', 'CDD:253677', 'CDD:276194', 'CDD:100107', 'CDD:100107'], ['GI:59800511', 'CDD:236307', 'CDD:249744', 'CDD:190425', 'CDD:214878'], ['GI:59800512', 'CDD:235777', 'CDD:133459', 'CDD:133459', 'CDD:133459'], ['GI:59800513'], ['GI:59800514', 'CDD:234666'], ['GI:59800515', 'CDD:235393', 'CDD:249744', 'CDD:190425', 'CDD:198164', 'CDD:249744', 'CDD:276202', 'CDD:238712', 'CDD:238712', 'CDD:238712', 'CDD:238712', 'CDD:238712'], ['GI:59800516'], ['GI:59800517', 'CDD:238509', 'CDD:238509'], ['GI:59800518', 'CDD:237139', 'CDD:198165', 'CDD:153215', 'CDD:153215', 'CDD:153215'], ['GI:59800520', 'CDD:225953', 'CDD:253258'], ['GI:59800521', 'CDD:183270', 'CDD:238994', 'CDD:238994'], ['GI:59800522', 'CDD:239309', 'CDD:239309'], ['GI:59800523', 'CDD:188209', 'CDD:197670', 'CDD:238042', 'CDD:238042', 'CDD:238042'], ['GI:59800524', 'CDD:131349'], ['GI:59800525', 'CDD:224129', 'CDD:133044', 'CDD:133044', 'CDD:133044'], ['GI:59800526', 'CDD:238266', 'CDD:238266', 'CDD:238266'], ['GI:59800527', 'CDD:236583', 'CDD:206687', 'CDD:206687', 'CDD:206687', 'CDD:206687', 'CDD:206687', 'CDD:206687', 'CDD:206687', 'CDD:206687', 'CDD:206687', 'CDD:133440'], ['GI:59800528'], ['GI:59800529', 'CDD:224748', 'CDD:260768', 'CDD:238141'], ['GI:59800530', 'CDD:235645', 'CDD:173902', 'CDD:173902', 'CDD:173902', 'CDD:173902', 'CDD:173902'], ['GI:651851637', 'CDD:276322'], ['GI:59800531', 'CDD:235536', 'CDD:185679', 'CDD:185679', 'CDD:214901'], ['GI:59800532', 'CDD:235588', 'CDD:275460', 'CDD:173912', 'CDD:173912', 'CDD:173912', 'CDD:173909', 'CDD:173909', 'CDD:173909', 'CDD:153414', 'CDD:153414', 'CDD:153414', 'CDD:253933'], ['GI:59800533', 'CDD:276322'], ['GI:59800534', 'CDD:234739'], ['GI:59800535', 'CDD:260116', 'CDD:260116', 'CDD:260116', 'CDD:260116'], ['GI:59800536', 'CDD:223620', 'CDD:119389', 'CDD:119389', 'CDD:257599'], ['GI:59800539', 'CDD:235554', 'CDD:213988', 'CDD:213988', 'CDD:213988', 'CDD:254393', 'CDD:258717', 'CDD:239931', 'CDD:239931', 'CDD:239931'], ['GI:59800540', 'CDD:238095', 'CDD:238095'], ['GI:59800541'], ['GI:59800542', 'CDD:223686', 'CDD:213993', 'CDD:213993'], ['GI:59800543', 'CDD:236517', 'CDD:99734', 'CDD:99734', 'CDD:99734', 'CDD:99734'], ['GI:59800544', 'CDD:224011', 'CDD:271430', 'CDD:187548', 'CDD:187548', 'CDD:187548', 'CDD:187548', 'CDD:187548'], ['GI:59800545', 'CDD:99740', 'CDD:99740', 'CDD:99740', 'CDD:99740'], ['GI:59800546', 'CDD:251270', 'CDD:100050', 'CDD:234265', 'CDD:100050', 'CDD:100050'], ['GI:59800547', 'CDD:275912', 'CDD:223515'], ['GI:59800548', 'CDD:223515', 'CDD:275912'], ['GI:59800549', 'CDD:225153', 'CDD:263758'], ['GI:59800550', 'CDD:232920', 'CDD:238611', 'CDD:238611', 'CDD:238611', 'CDD:224896'], ['GI:59800551', 'CDD:234774', 'CDD:251987'], ['GI:59800552', 'CDD:238889', 'CDD:238889', 'CDD:249761', 'CDD:238889'], ['GI:59800553', 'CDD:173954', 'CDD:173954', 'CDD:173954', 'CDD:173954'], ['GI:59800554', 'CDD:223775', 'CDD:238260', 'CDD:238260', 'CDD:238260', 'CDD:238260'], ['GI:59800555', 'CDD:256588', 'CDD:233905', 'CDD:198033', 'CDD:252276', 'CDD:249725'], ['GI:59800556', 'CDD:225709'], ['GI:59800557', 'CDD:260846'], ['GI:59800558', 'CDD:225707'], ['GI:59800559', 'CDD:227306'], ['GI:59800560', 'CDD:227342', 'CDD:250219', 'CDD:250215'], ['GI:59800561', 'CDD:223296', 'CDD:206665', 'CDD:206665', 'CDD:206665', 'CDD:206665', 'CDD:206665', 'CDD:206665', 'CDD:206665', 'CDD:206665', 'CDD:206665'], ['GI:59800562', 'CDD:225418', 'CDD:249527'], ['GI:59800563', 'CDD:224252'], ['GI:59800564', 'CDD:275607'], ['GI:59800565', 'CDD:236585', 'CDD:276486'], ['GI:59800566', 'CDD:153246', 'CDD:153246'], ['GI:161572979', 'CDD:235374', 'CDD:249858', 'CDD:251520', 'CDD:111646'], ['GI:59800568', 'CDD:224938', 'CDD:251093'], ['GI:59800570'], ['GI:59800571', 'CDD:225407', 'CDD:252049'], ['GI:59800572', 'CDD:223715', 'CDD:119399', 'CDD:119399', 'CDD:119399', 'CDD:238030', 'CDD:238030', 'CDD:238030', 'CDD:238030'], ['GI:59800573', 'CDD:183285', 'CDD:239900', 'CDD:239900', 'CDD:239900', 'CDD:255792'], ['GI:59800574', 'CDD:239510', 'CDD:239510', 'CDD:239510'], ['GI:59800576', 'CDD:213561', 'CDD:238312', 'CDD:238312'], ['GI:59800577', 'CDD:236794', 'CDD:239934', 'CDD:239934', 'CDD:239934', 'CDD:238005', 'CDD:238005', 'CDD:238005', 'CDD:238034', 'CDD:238034', 'CDD:238034'], ['GI:59800578', 'CDD:234761', 'CDD:250373', 'CDD:187535', 'CDD:268288'], ['GI:59800579', 'CDD:179791'], ['GI:59800580'], ['GI:59800581', 'CDD:225406'], ['GI:59800582'], ['GI:59800583', 'CDD:130902', 'CDD:271629', 'CDD:271712', 'CDD:99778', 'CDD:99778', 'CDD:99778', 'CDD:99778', 'CDD:271757', 'CDD:257762'], ['GI:59800584'], ['GI:59800585'], ['GI:745998704', 'CDD:239100', 'CDD:239100', 'CDD:226359', 'CDD:276299', 'CDD:213179', 'CDD:213179', 'CDD:213179'], ['GI:59800586', 'CDD:179791', 'CDD:273260'], ['GI:59800588', 'CDD:238953', 'CDD:238953'], ['GI:59800589', 'CDD:148918'], ['GI:59800590', 'CDD:223792', 'CDD:250470', 'CDD:238351'], ['GI:59800591', 'CDD:226282', 'CDD:223675', 'CDD:100105', 'CDD:100105'], ['GI:59800593', 'CDD:235417'], ['GI:59800594', 'CDD:224671', 'CDD:252025'], ['GI:59800595'], ['GI:59800596', 'CDD:233695', 'CDD:257693', 'CDD:238487', 'CDD:238487', 'CDD:238487', 'CDD:238487'], ['GI:59800597', 'CDD:182661', 'CDD:238013', 'CDD:238013', 'CDD:238013', 'CDD:238013', 'CDD:238013', 'CDD:197771'], ['GI:59800598', 'CDD:224159'], ['GI:59800599'], ['GI:59800600', 'CDD:223809', 'CDD:119392', 'CDD:119392'], ['GI:59800601', 'CDD:224671', 'CDD:252025'], ['GI:59800602', 'CDD:226910', 'CDD:239963', 'CDD:215020'], ['GI:59800603', 'CDD:223396'], ['GI:59800604', 'CDD:234612', 'CDD:270364', 'CDD:270364', 'CDD:270364'], ['GI:651851638', 'CDD:257686'], ['GI:59800607', 'CDD:238167', 'CDD:214692', 'CDD:238167', 'CDD:238167', 'CDD:238167', 'CDD:238034', 'CDD:238034', 'CDD:238034'], ['GI:59800609', 'CDD:129820', 'CDD:239200', 'CDD:239200', 'CDD:239200', 'CDD:239200', 'CDD:239200'], ['GI:651851625', 'CDD:179348'], ['GI:161572978', 'CDD:238294', 'CDD:238294', 'CDD:238294', 'CDD:238294'], ['GI:59800612', 'CDD:153246', 'CDD:153246'], ['GI:59800614', 'CDD:99825', 'CDD:99825', 'CDD:99825', 'CDD:99825', 'CDD:255760'], ['GI:59800615', 'CDD:237585', 'CDD:189007', 'CDD:189007', 'CDD:189007'], ['GI:59800616'], ['GI:59800617'], ['GI:59800619'], ['GI:59800620', 'CDD:271796'], ['GI:59800621'], ['GI:59800622'], ['GI:59800624', 'CDD:223874', 'CDD:275504', 'CDD:238347'], ['GI:59800625', 'CDD:224033', 'CDD:275585', 'CDD:119348', 'CDD:119348', 'CDD:119348'], ['GI:59800626', 'CDD:224046', 'CDD:213202'], ['GI:59800627', 'CDD:235418'], ['GI:59800628', 'CDD:234581'], ['GI:59800629', 'CDD:223877', 'CDD:250860', 'CDD:253112'], ['GI:59800630', 'CDD:234588'], ['GI:59800631', 'CDD:223715', 'CDD:100122', 'CDD:100122', 'CDD:119399', 'CDD:119399', 'CDD:119399', 'CDD:238030', 'CDD:238030', 'CDD:238030', 'CDD:238030'], ['GI:59800632', 'CDD:223816', 'CDD:238088', 'CDD:238088', 'CDD:238088', 'CDD:238088', 'CDD:238088', 'CDD:238225', 'CDD:238225'], ['GI:59800633', 'CDD:225844', 'CDD:252887', 'CDD:256667'], ['GI:59800634', 'CDD:242419'], ['GI:59800635', 'CDD:238310', 'CDD:238310', 'CDD:238310'], ['GI:59800636', 'CDD:275617'], ['GI:59800637', 'CDD:234818', 'CDD:275655'], ['GI:59800638', 'CDD:235091'], ['GI:59800639', 'CDD:238607', 'CDD:238607', 'CDD:238607', 'CDD:238607'], ['GI:59800640', 'CDD:212141', 'CDD:212141'], ['GI:59800641', 'CDD:223991', 'CDD:176195', 'CDD:176195', 'CDD:176195'], ['GI:59800642', 'CDD:212541', 'CDD:212541', 'CDD:212541'], ['GI:59800643', 'CDD:251482'], ['GI:59800644', 'CDD:258018', 'CDD:235615', 'CDD:254281', 'CDD:273336'], ['GI:59800645', 'CDD:237275', 'CDD:254281', 'CDD:216990'], ['GI:59800646', 'CDD:238213', 'CDD:238213', 'CDD:238213', 'CDD:238213'], ['GI:59800647', 'CDD:183226', 'CDD:276299', 'CDD:213179', 'CDD:213179', 'CDD:213179', 'CDD:213179', 'CDD:213179', 'CDD:213179', 'CDD:213179', 'CDD:254779'], ['GI:59800648', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394'], ['GI:59800649', 'CDD:224098', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394'], ['GI:59800650', 'CDD:276298', 'CDD:223737'], ['GI:59800651', 'CDD:223083'], ['GI:59800652', 'CDD:236490', 'CDD:198027', 'CDD:239906', 'CDD:239906', 'CDD:238548', 'CDD:238548', 'CDD:238548', 'CDD:238548', 'CDD:238548'], ['GI:59800653', 'CDD:235809', 'CDD:250534', 'CDD:249825', 'CDD:271901'], ['GI:59800654', 'CDD:226201'], ['GI:59800655', 'CDD:224719'], ['GI:59800656', 'CDD:276303', 'CDD:223620', 'CDD:119389', 'CDD:119389'], ['GI:59800657', 'CDD:153097', 'CDD:153097'], ['GI:59800658', 'CDD:178807'], ['GI:59800659', 'CDD:270377', 'CDD:270377', 'CDD:270377'], ['GI:59800660', 'CDD:223539', 'CDD:132997', 'CDD:132997'], ['GI:59800661', 'CDD:223735', 'CDD:100051', 'CDD:100051', 'CDD:100051'], ['GI:59800662', 'CDD:179105', 'CDD:206732', 'CDD:206732', 'CDD:206732', 'CDD:206732', 'CDD:206732', 'CDD:206732', 'CDD:206732', 'CDD:206732', 'CDD:206732', 'CDD:206732', 'CDD:239660', 'CDD:258632'], ['GI:59800663', 'CDD:234598'], ['GI:59800664', 'CDD:240082', 'CDD:240082', 'CDD:240082'], ['GI:59800665', 'CDD:240083', 'CDD:240083'], ['GI:59800666', 'CDD:153219', 'CDD:171871', 'CDD:153219', 'CDD:153219', 'CDD:153219'], ['GI:59800667', 'CDD:223357'], ['GI:59800668', 'CDD:226361', 'CDD:213226', 'CDD:213226', 'CDD:213226', 'CDD:213226', 'CDD:213226', 'CDD:213226', 'CDD:213226', 'CDD:213226'], ['GI:59800669', 'CDD:224099', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394', 'CDD:119394'], ['GI:59800670', 'CDD:270261', 'CDD:270261'], ['GI:59800672', 'CDD:179143', 'CDD:176463', 'CDD:176463', 'CDD:176463'], ['GI:59800673', 'CDD:133021', 'CDD:182252', 'CDD:133021', 'CDD:133021'], ['GI:59800674', 'CDD:238285', 'CDD:238285', 'CDD:238285'], ['GI:59800675'], ['GI:59800676', 'CDD:238239', 'CDD:238239', 'CDD:238239', 'CDD:238239'], ['GI:59800677', 'CDD:240022', 'CDD:240022'], ['GI:59800678', 'CDD:191481', 'CDD:258590'], ['GI:59800679'], ['GI:59800680', 'CDD:257184'], ['GI:59800681'], ['GI:59800682'], ['GI:59800683', 'CDD:223246', 'CDD:266653'], ['GI:59800684', 'CDD:234673', 'CDD:163665', 'CDD:163665', 'CDD:163665'], ['GI:59800685', 'CDD:100004', 'CDD:100004', 'CDD:100004'], ['GI:59800686', 'CDD:111368'], ['GI:59800687', 'CDD:236346', 'CDD:100105', 'CDD:100105', 'CDD:253997'], ['GI:59800688', 'CDD:236137', 'CDD:238062', 'CDD:238062', 'CDD:238062', 'CDD:238062', 'CDD:145978', 'CDD:251739', 'CDD:257328', 'CDD:237994', 'CDD:237994', 'CDD:237994'], ['GI:59800689', 'CDD:225629', 'CDD:129010'], ['GI:59800690', 'CDD:225567', 'CDD:133475', 'CDD:133475', 'CDD:133475', 'CDD:133475'], ['GI:59800692', 'CDD:234814', 'CDD:238835', 'CDD:238835', 'CDD:238835'], ['GI:59800693', 'CDD:223358', 'CDD:249824', 'CDD:133453', 'CDD:133453'], ['GI:59800695', 'CDD:234808']]
            print(tabela_CDD(GI))
        elif ans=="60":
            ans = False
        else:
            print("\nInvalid")



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
    


#main
if __name__ == "__main__":
    res=str(input(print("Ja tem ficheiro S/N? "))).upper()
    if res=="N":
        record=create_file()
    elif res=="S":
        filename=str(input("Nome do ficheiro? "))+".gb"
        record = SeqIO.read("../res/"+filename, "genbank") 
    menu(record)
    
