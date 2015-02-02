class Record(dict):
    """
    This record stores the information of one keyword or category in the
    keywlist.txt as a Python dictionary. The keys in this dictionary are
    the line codes that can appear in the keywlist.txt file:

    ---------  ---------------------------     ----------------------
    Line code  Content                         Occurrence in an entry
    ---------  ---------------------------     ----------------------
    ID         Identifier (keyword)            Once; starts a keyword entry.
    AC         Accession (KW-xxxx)             Once.
    DE         Definition                      Once or more.
    CC         Subcellular Location            Once or more; comments.
    SQ         Sequence                        Once; contains only the heading information.
    """
    def __init__(self):
        dict.__init__(self)
        for keyword in ("DE", "CC"):
            self[keyword] = []
    
def parse(handle): # The parameter handle is the UniProt KeyList file.
    record = Record()
    # Now parse the records
    for line in handle:
        key = line[:2]
        if key=="//": # The last line of the current record has been reached.
            record["DE"] = " ".join(record["DE"])
            record["CC"] = " ".join(record["CC"])
            yield record     # So we output the record and pass to other one. 
            record = Record()
        elif line[2:5]=="   ": # If not, we continue recruiting the information. 
            value = line[5:].strip()
            if key in ("ID", "AC", "SQ"):
                record[key] = value
            elif key in ("DE", "CC"):
                record[key].append(value) 
            else:
                pass
            
    # Read the footer and throw it away
    for line in handle:
        pass