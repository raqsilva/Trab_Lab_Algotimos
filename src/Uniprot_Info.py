# -*- coding: utf-8 -*-
import urllib #getting info from site
from Uniprot_Parser import * #parsing uniprot text file
from Bio.SeqIO import UniprotIO #parsing uniprot xml file


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
    