import matplotlib.pyplot as plt
import cPickle as pickle

def peptide_set(s, k, l, offset=0):
    if (k==1):
        return s
    if (len(s) < k):
        return []
    i = offset
    plist = []
    while (i+k <= len(s)):
        plist.append(s[i:i+k])
        i = i + l
    return plist

def label_set(s, k, l, offset=0):
    if (k==1):
        return s
    if (len(s) < k):
        return []
    i = offset
    plist = []
    while (i+k <= len(s)):
        plist.append(s[i:i+k])
        i = i + l
    return plist


def show_data(master_records, K=15, L=1, source="AGGREGATE"):
    num_epitopes = 0
    X=[]
    y=[]
    # model each amino acid position using a 15-mer, set up training examples as positive/negative
    # by whether central AA is contained in an epitope
    for p in master_records.keys():
        #print p
        # exclude this one for testing purposes
        antigen_seq = master_records[p][0]
        num_epitopes += len(master_records[p][1])
        curr_epitopes = master_records[p][1]
        #print curr_epitopes
        # WE NEED TO DEAL WITH THE END OF THE SEQUENCE... WE MIGHT HAVE A LOT OF MISSING POSITIVES....
        curr_peptides = peptide_set(antigen_seq, K, L)
        # 0: sequence conservation; 1: b-factor; 2: COREX; 3: solvent accessible area; 4: aggregate
        if (source == "SEQUENCE"):
            curr_stability_data = peptide_set(master_records[p][2][0], K, L)
        elif (source == "BFACTOR"):
            curr_stability_data = peptide_set(master_records[p][2][1], K, L)
        elif (source == "COREX"):
            curr_stability_data = peptide_set(master_records[p][2][2], K, L)
        elif (source == "ASA"):
            curr_stability_data = peptide_set(master_records[p][2][3], K, L)
        else: #source should be "AGGREGATE"
            curr_stability_data = peptide_set(master_records[p][2][4], K, L)
        # get all stability data individually
        #curr_stability_data = zip(peptide_set(master_records[p][2][0], K, L), peptide_set(master_records[p][2][1], K, L), peptide_set(master_records[p][2][2], K, L), peptide_set(master_records[p][2][3], K, L))
        #curr_epitope_mask = peptide_set(master_records[p][3], K, L)
        curr_epitope_mask=(label_set(master_records[p][3],K,L))

        for i in range(len(curr_stability_data)):
            X.append(curr_stability_data[i])
            y.append(curr_epitope_mask[i])

    print "We get a total of %d %d-mers"%(len(X),K)
    print "\nEach %d-mer is a training instance"%(K)


master_records = pickle.load(open('3_30_human.pickle', 'r'))
print "loaded %d protein data sets as follows\n" % len(master_records)

num_epitopes = sum([len(x[2]) for x in master_records.values()])

print master_records.keys()

print "\ntotal number of eptiopes is %d\n" %(num_epitopes)


show_data(master_records, K=15, L=1, source="AGGREGATE")






