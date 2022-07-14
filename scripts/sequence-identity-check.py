import numpy as np
import os
import sys
import csv
import glob
from operator import itemgetter

import Bio
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist

#Jonathan Borowsky
#7/14/22
#this program calculates the sequence identity between each pair of proteins
#that passed manual filtering and returns a new protein list without the
#lower-sequence-identity member of each pair of overlapping proteins removed

################################################################################

directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/cryptic-pocket-filtering-3"

#keep track of different versions of input and output files
#(e.g. one for forward and one for reverse pockets)
serial_in = 1
serial_out = 1

#consider checking sequence identity before filtering manually
manually_filtered_proteins = np.load(f"{directory}/iofiles-manual/screened-pairs-v{serial_in}.npy")

#add an option to download fasta files for the apo structures, add a list of pdb ids, or draw from a directory of fasta files

proteins_all = [e[0] for e in manually_filtered_proteins]
print(f"Holo PDB IDs: {proteins_all}")

################################################################################
#check sequence identity of each pair of proteins
################################################################################


id_cutoff = 0.4 #maximum sequence identity

#define which RMSD to use when automatically removing one protein from a pair
# with sequence identity exceeding the cutoff
rmsd_types_inds = {"all c alpha": -4, "cryptic site": -3}
rmsd_type = "all c alpha"
rmsd_ind = rmsd_types_inds[rmsd_type]

nprots = len(proteins_all)

zmax = nprots*(nprots-1)/2
z = 0

#store the sequence identity of every protein pair
seqid_matrix = np.zeros([nprots, nprots])

#collect a list of matches above 40% identity
seqid_out = []
x_removal = [] #indices of structures to remove

print(f"aligning {len(proteins_all)} proteins; {int(zmax)} protein pairs")

for x, prot1 in enumerate(proteins_all):

    print(f"{x}: {prot1}, {'{:.4f}'.format(z/zmax)} complete")

    fastaname1 = f"{directory}/iofiles-blast/rcsb-fasta-ab/{prot1.upper()}.fasta"
    try:
        fasta_seq1 = SeqIO.read(open(fastaname1),'fasta').seq
    except:
        print("issue reading sequence, please provide a cleaned fasta file:")
        print(SeqIO.read(open(fastaname1),'fasta'))

    for y, prot2 in enumerate(proteins_all[x+1:]):

        fastaname2 = f"{directory}/iofiles-blast/rcsb-fasta-ab/{prot2.upper()}.fasta"
        try:
            fasta_seq2 = SeqIO.read(open(fastaname2),'fasta').seq
        except:
            print("issue reading sequence, please provide a cleaned fasta file:")
            print(SeqIO.read(open(fastaname2), 'fasta'))


        alignment = pairwise2.align.globalds(fasta_seq1, fasta_seq2, matlist.blosum62, -11, -1) #parameters matching blastp
        #setting d uses the matrix, s uses the gap creation and extension penalties
        #note that format_alignment() can be used to display alignments nicely as long as they aren't wider than the terminal

        #print(alignment)

        minlen = min(len(fasta_seq1), len(fasta_seq2))
        #use the shorter of the two sequences as the denominator to determine the percent identity
        #as a fraction of the maximum achieveable percent identity given the sequence lengths
        #this approach is more conservative

        identical_resis = 0

        for k in range(len(alignment[0].seqA)):
            if alignment[0].seqA[k] == alignment[0].seqB[k]:
                identical_resis += 1

        #qc, remove this
        if len(alignment[0].seqA)/len(alignment[0].seqB) != 1:
            print("----------------------------------------------------------------\n---------------------------------------eep------------------------------")

        seqid = identical_resis/minlen #fraction of identical residues in aligned sequences

        seqid_matrix[x][y] = seqid

        #only works well if the sequence space is only sparsely populated, which it appears to be for the set of filtered pairs
        if seqid > id_cutoff:

            seqid_out.append([seqid, prot1, prot2])

            rmsd_x = manually_filtered_proteins[x][rmsd_ind]
            rmsd_y = manually_filtered_proteins[y][rmsd_ind]

            if rmsd_x < rmsd_y:
                print(f"removing {manually_filtered_proteins[x]}")
                x_removal.append(x)
            elif rmsd_x > rmsd_y:
                print(f"removing {manually_filtered_proteins[y]}")
                x_removal.append(y)
            elif rmsd_x == rmsd_y:
                print(f"Proteins have equal RMSD; removing {manually_filtered_proteins[y]} arbitrarily.")
                x_removal.append(y)

            print(f"{prot1} and {prot2} have {seqid*100}% identity")

        z+=1 #counter

#save PDB ID list and sequence identity matrix
np.save(f"{directory}/iofiles-seqid/holo-pdb-ids-v{serial_out}", proteins_all)
np.save(f"{directory}/iofiles-seqid/sequence-identity-matrix-v{serial_out}", seqid_matrix)

seqid_out_sorted = sorted(seqid_out, key=itemgetter(0), reverse = True)
np.save(f"{directory}/iofiles-seqid/sequence-id-above-{id_cutoff}-v{serial_out}", seqid_out_sorted)

#-------------------------------------------------------------------------------
#save an updated npy list and csv of remaining proteins

screened_pairs_seqid = [e for x, e in enumerate(manually_filtered_proteins) if x not in x_removal]
np.save(f"{directory}/iofiles-seqid/screened-pairs-seqid-v{serial_out}", screened_pairs_seqid)

# open the file in write mode
f = open(f"{directory}/iofiles-seqid/screened-pairs-seqid-v{serial_out}", 'w', newline='')
# create the csv writer
writer = csv.writer(f, quoting = csv.QUOTE_NONE)

for pair in screened_pairs_seqid:
    print(pair)
    writer.writerow(pair)

#-------------------------------------------------------------------------------

#print list of protein pairs with excessive sequence identity
print("apo-holo pairs with sequence identities above 40%:")
print(seqid_out_sorted)

#print(f"sequence identity matrix: {seqid_matrix}")
