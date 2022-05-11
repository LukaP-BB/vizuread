import sys

"""
Script parsant le fichier all_trans_no_filter et affichant les translocs à vérifier
""" 

if len(sys.argv) != 2 :
    print("Usage : python3 parse_alltrans.py /path/to/all_trans.txt")
    exit(1)

with open(sys.argv[1]) as fp :
    
    dico = {}
    
    for line in fp :
        fields = line.split("\t")
        if len(fields) == 1 :
            current_id = fields[0].strip()
            dico[current_id] = {}
        
        if len(fields) == 9 :
            current_transloc = fields[0]
            dico[current_id][current_transloc] = []
            
            reads, total_reads = int(fields[-2]), int(fields[-1])
            dico[current_id][current_transloc].append((reads, total_reads, fields[1:4]))

        if len(fields) == 8 :
            reads, total_reads = int(fields[-2]), int(fields[-1])
            dico[current_id][current_transloc].append((reads, total_reads, fields[0:3]))



for patient, translocs in dico.items() :
    to_print = False
    to_verify = []
    for transloc, read_infos in translocs.items() :
        is_already_accepted = False
        for infos in read_infos :
            reads, total_reads, pos = infos
            if reads >= 5 and total_reads >= 10 :
                is_already_accepted = True
                break
            if total_reads < 10 : 
                to_print = True
                to_verify.append(f"{transloc.ljust(8)}: {' '.join(pos).rjust(25)} | {reads, total_reads}")

    if to_print :
        print(patient + " : \n\t" + "\n\t".join(to_verify))