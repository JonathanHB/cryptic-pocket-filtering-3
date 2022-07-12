import re
import numpy as np
import os
import Bio.Blast
from Bio.Blast import NCBIXML

#take the lists of PDB structures downloaded from MOAD (https://www.bindingmoad.org/Home/download)
#for each structure:
    #download the fasta sequence
    #blastp it against PDB structures of at least 50 nucleotides
    #save the resulting structures with single hsps 100% sequence identity and sufficient coverage
    #   note that having a coverage requirement below 1 leaves open the possibility of mutations outside the hsp
    #   but having a coverage requirement of 1 excludes many identical structures of slightly different sequence lengths
    #   (i.e. a disordered terminal tail might be truncated in different places in different structures)

#path to filtering directory
directory = "/project/bowmanlab/borowsky.jonathan/FAST-cs/protein-sets/cryptic-pocket-filtering-3" #new_pockets/iofiles"
#part = "all" #part a or b or both ("all") of the MOAD database, which is broken in two on account of its large size

def makedirectories():
    #os.system(f"cd {directory}/iofiles-blast")
    subdirectories = ["rcsb-fasta-ab", "blast-xml-out", "blast-output", "blast-output-ho", "blast-output-m", "blast-output-m-ho"]
    #"moad_xml", "rcsb_pdb_holo", "rcsb_pdb_apo", "monomer_apo", "monomer_holo", "pairing_index", "note_files", "qc_structures"]
    for sdir in subdirectories:
        os.system(f"mkdir {directory}/iofiles-blast/{sdir}")

makedirectories()

coverage_limit = 0.90 #minimum fractional overlap of allowable hits; filter_blast2pairs.py subsequently makes sure no internal mismatches get through this way

#get the PDB IDs of all the MOAD structures
structures = np.hstack([np.unique([i[0:4].upper() for i in os.listdir(f"{directory}/iofiles-blast/every_part_a/BindingMOAD_2020/")]), np.unique([i[0:4].upper() for i in os.listdir(f"{directory}/iofiles-blast/every_part_b/BindingMOAD_2020/")])])
n_structures = len(structures)

initial_ind = 0
final_ind = 200#n_structures
print(f"processing {final_ind-initial_ind} structures from indices {initial_ind} to {final_ind}.")

#structures for which fasta files have already been downloaded
existing_ids = [i[0:4] for i in os.listdir(f"{directory}/iofiles-blast/rcsb-fasta-ab/")]

#loop through moad structures
for i in structures[initial_ind:final_ind]:

    hits = []
    print(i)

    #download fasta files if necessary
    if i not in existing_ids:
        os.system(f"wget https://www.rcsb.org/fasta/entry/{i}/display -O {directory}/iofiles-blast/rcsb-fasta-ab/{i}.fasta")

    #blast the MOAD sequence against PDB
    #note that the ftp version of blastp must be used because the apt-get version is out of date and can't load the pdbaa database.
    #-outfmt 5 is required to produce a biopython-readable .xml file
    os.system(f"/project/bowmanlab/borowsky.jonathan/installations/ncbi-blast-2.12.0+/bin/blastp -db /project/bowmanlab/borowsky.jonathan/installations/pdbaa -query {directory}/iofiles-blast/rcsb-fasta-ab/{i}.fasta -outfmt 5 -out {directory}/iofiles-blast/blast-xml-out/{i}-results.xml")

    #load the blastp results
    result_handle = open(f"{directory}/iofiles-blast/blast-xml-out/{i}-results.xml")

    #skip structures containing multiple different protein sequences
    #As multiple proteins are usually added to crystallization solutions for the purpose of studying complexes thereof,
    #there are likely relatively few multi protein PDB structures with monomeric bioassemblies
    #so I'm skipping them for the time being for ease of processing
    #The following can be used to read xml results for multimers, but further downstream adjustments would be required
    #  "blast_records = NCBIXML.parse(result_handle)
    #   blast_record = next(blast_records)"
    try:
        blast_record = NCBIXML.read(result_handle)
    except ValueError as err:
        print("Encountered an error, likely due to the presence of multiple different sequences in the fasta file.\
        This structure will be skipped because this filtering pipeline is designed to process monomeric structures.\
        The error is returned below:")
        print(err)
        continue

    #for each structure aligned to the moad structure
    for alignment in blast_record.alignments:
        print(alignment.title)
        hsp = alignment.hsps[0] #there should be only one hsp for any alignment good enough to be useful; this is checked below
        #print(hsp.identities)
        #check that there is one hsp with 100% identity and good coverage
        if len(alignment.hsps) == 1 and hsp.identities == hsp.align_length and hsp.align_length/alignment.length > coverage_limit:

            entries = re.split("pdb\||>pdb\|", alignment.title) #separate out all pdb files with a given sequence
            pdb_ids = [title[0:4] for title in entries[1:]] #the first entry is "" from the left of the first pdb| header
            hits = hits+pdb_ids

    #make the query id, if any, the first element in the array (I think it is before np.unique() rearranges things,
    #but I'm not sure and this way if the blast results are in an unusual order it still works)

    hits_u = np.unique(hits)

    hits_apo_a  = [] #all non-moad hits
    hits_apo_m = np.delete(hits_u, np.where(i.upper()==hits_u)) #all hits except the holo structure
    for acand in hits_u: #for all unique hits
        if acand not in structures:
            hits_apo_a.append(acand) #add all non-moad structures

    #add the holo structure id to the start of the list
    hits_a = np.hstack([i.upper(), hits_apo_a])
    hits_m = np.hstack([i.upper(), hits_apo_m])

    print(hits_a)

    #save the output with and without structures with no apo candidates
    if len(hits_a) > 1:
        np.save(f"{directory}/iofiles-blast/blast-output/{i}-hits.npy", hits_a)
    else:
        np.save(f"{directory}/iofiles-blast/blast-output-ho/{i}-hits-ho.npy", hits_a)

    if len(hits_m) > 1:
        np.save(f"{directory}/iofiles-blast/blast-output-m/{i}-hits-m.npy", hits_m)
    else:
        np.save(f"{directory}/iofiles-blast/blast-output-m-ho/{i}-hits-m-ho.npy", hits_m)

#---------------------------------------------------End of Code---------------------------------------------------

#this script ran successfully from 8/15/21 to 8/16/21, processing all of MOAD
