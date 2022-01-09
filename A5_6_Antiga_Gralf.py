#The assignment was discussed with my colleagues Oliver Strauss and Thresa Plumer
#https://github.com/Johlf/A5_6_Antiga_Gralf
import Bio
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO, Seq, Align, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.SeqUtils import GC
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.SeqUtils.ProtParam import ProteinAnalysis
Entrez.email = 'be19b054@technikum-wien.at'

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Strg+F8 to toggle the breakpoint.

if __name__ == '__main__':
    # https://biopython.org/docs/1.76/api/Bio.Entrez.html
    search_string = "SNCA"
    handle = Entrez.esearch(db="nucleotide", term=search_string, retmax=69)
    record = Entrez.read(handle)
    ids = record['IdList']
    handle.close()
    print(ids)

    # goes through every id and checks if the species is already saved in orthologs and saves new species in
    # the orthologs and orthorec list.
    orthologs = []
    orthrec = []
    for x in ids:
        handle = Entrez.efetch(db="nucleotide", id=x, rettype="gb", retnode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        # Oliver Strauss inspired me a bit in this section.
        temp = record.annotations["organism"]
        if not temp in orthologs:
            orthrec.append(record)
            orthologs.append(temp)
            print(temp)

        if len(orthologs) == 5:
            break

    print(orthrec)
    f = open("standard.csv", "w+")
    data = pd.read_csv('standard.csv', sep=";", header=None,names=["Accession Number", "Title", "Organism", "Sequence length", "Sequence"])
    df=pd.DataFrame(data)
    accessionnr = []
    title = []
    organism = []
    length = []
    seq = []
    for x in range (len(orthrec)):
        accessionnr.append(orthrec[x].id)
        df.loc[x, 'Accession Number'] = orthrec[x].id
        title.append(orthrec[x].description)
        df.loc[x, 'Title'] = orthrec[x].description
        organism.append(orthrec[x].annotations["organism"])
        df.loc[x, 'Organism'] = orthrec[x].annotations["organism"]
        length.append(len(orthrec[x].seq))
        df.loc[x, 'Sequence length'] = str(len(orthrec[x].seq))
        seq.append(orthrec[x].seq)
        df.loc[x, 'Sequence'] = orthrec[x].seq
    print(df)

    df.to_csv("standard.csv", sep = ";")
    #https: // biopython.org / docs / 1.75 / api / Bio.SeqUtils.html
    GCPercentage=[]
    for x in range(len(seq)):
        GCPercentage.append(GC(seq[x]))
    print(GCPercentage)
    #https://stackoverflow.com/questions/32833230/biopython-alignio-valueerror-says-strings-must-be-same-length
    #makes sequence same length, otherwise error during alignment
    maxlen = max(len(record.seq) for record in orthrec)
    for record in orthrec:
        if len(record.seq) != maxlen:
            sequence = str(record.seq).ljust(maxlen, '.')
            record.seq = Seq.Seq(sequence)
    assert all(len(record.seq) == maxlen for record in orthrec)

    #https://biopython.org/docs/1.76/api/Bio.Align.html#subpackages
    align = MultipleSeqAlignment(orthrec)
    print(align)

    #https: // www.youtube.com / watch?v = wBdz3vFQ4Ks
    #calculates the distances between the sequences
    calculator = DistanceCalculator("identity")
    distance_matrix=calculator.get_distance(align)
    print(distance_matrix)
    # generates a phylogenetic tree
    constructor=DistanceTreeConstructor(calculator)
    tree=constructor.build_tree((align))
    tree.rooted=True
    print(tree)
    Phylo.write(tree,"tree.xml","phyloxml")
    fig=Phylo.draw(tree)

    # translate RNA to Protein
    Protein=[]
    for x in range(len(seq)):
        Protein.append(seq[x].translate(to_stop=True))
    print(Protein)

    # https: // biopython.org / docs / 1.76 / api / Bio.SeqUtils.ProtParam.html
    # calculates instability, aromaticity, isoelectric point and analyses the structure
    Instability=[]
    Aromaticity=[]
    Isoelectric_point=[]
    Secstruc=[]
    for x in range(len(Protein)):
        X=ProteinAnalysis(str(Protein[x]))
        Instability.append("%0.2f" % X.instability_index())
        Aromaticity.append("%0.2f" % X.aromaticity())
        Isoelectric_point.append("%0.2f" % X.isoelectric_point())
        sec_struc=X.secondary_structure_fraction()
        Secstruc.append(sec_struc)

    #adds instability and aromaticity to the already existing table
    df['Instability']= Instability
    df['Aromaticity'] = Aromaticity
    df.to_csv("standard.csv", sep =";")
    print(df)