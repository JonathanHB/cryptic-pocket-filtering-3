import numpy as np
import os

#take the pairs saved to pairing-output/ and sort them in descending order by minimum rmsd to each holo structure
#this produces a list which begins with pairs with the largest conformational change not captured by experimental apo variation
#when multiple holo structures have the same apo structure, only the highest-rmsd pair is returned
#experimental holo variation can be attributable to differences in ligand size, and this code does not attempt to separate this from other factors

manual_removal = []
#original value of manual removal and its description:
#["5x93", "1q6m", "1q6j", "4acx", "4b70", "4b05", "1yp9", "3ufl"]
#structures that are messed up in rare ways such that they fail silently and have to be removed manually
#5x93 has an incomplete remark 465
#"1q6m", "1q6j", "4acx", 4b70", 1yp9, 3ufl, and 4b05 have their residues numbered out of order
#in the case of 4b05 this also causes the gap-checker to fail since it believes
#that the remark 465 residues are all before the start of the structure after reading 498 as the first residue
#2wkp and 1opk are each missing half the structure

directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/cryptic-pocket-filtering-3"

#-------------------------------------------------------------------------------
#compile all the holo-structure-level individual output files from filtering into one list

pair_ligands = []

directory1 = f"{directory}/iofiles-pdb"
holo_list1 = os.listdir(f"{directory1}/pairing-output")

for i in holo_list1:
    if i[0:6] == "ligand":
        pair = np.load(f"{directory1}/pairing-output/{i}", allow_pickle=True)
        if len(pair[0][5])!=0 and (pair[0][0] not in manual_removal):
            pair_ligands.append(pair)

#use this to merge filtering output files fron two different folders into a single list
mergeresults = False

if mergeresults:
    directory2 = f"{directory}/iofiles-pdb"
    holo_list2 = os.listdir(f"{directory2}/pairing-output")

    for i in holo_list2:
        if i[0:6] == "ligand":
            pair = np.load(f"{directory2}/pairing-output/{i}", allow_pickle=True)
            if len(pair[0][5])!=0 and (pair[0][0] not in manual_removal):
                pair_ligands.append(pair)

pair_ligands = sorted(pair_ligands, key = lambda x: x[0][4], reverse = True)

#-------------------------------------------------------------------------------
#copy the first (and therefore highest RMSD since the array is sorted)
#apo-holo pair with each apo or holo structure into a new array

apo_ids = []
holo_ids = []
lowrmsd_pairs = []
for j in pair_ligands:
    #avoid returning multiple holo structures with the same apo or holo structure
    #the case of duplicate holo structures should only occur if multiple sets of
    #results are being merged
    if j[0][2] not in apo_ids and j[0][0] not in holo_ids:
        lowrmsd_pairs.append(j[0])
        apo_ids.append(j[0][2])
        holo_ids.append(j[0][0])

#-------------------------------------------------------------------------------
#save output

serial = "1" #for tracking different versions
np.save(f"{directory}/filtering-output/all-apoholo-pairs-by-all-c-alpha-rmsd-v{serial}", np.array(lowrmsd_pairs, dtype = object))
