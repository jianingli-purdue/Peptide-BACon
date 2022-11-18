from schrodinger import structure
from schrodinger.structutils import analyze, build, transform, minimize
from schrodinger.structutils.build import add_hydrogens, delete_hydrogens
from math import pi
# from schrodinger.protein import buildpeptide
from schrodinger.protein import helix
import os
from schrodinger import structure
from schrodinger.infra import mm
from schrodinger.structutils import analyze, build, transform, minimize

##------------------------------------------------------------------------------------------------------------
# function to add correct terminals to alpha helix
# FIXME: need to change atomic weight, etc for new atoms

def fix_helix(st, seq):
    deleteme = []
    delete_hydrogens(st)
    for res in st.residue:
        if res.resnum == 1:
            for atom in res.atom:
                deleteme.append(atom.index)
        if res.resnum == 2:
            for atom in res.atom:
                if atom.pdbname.strip() == "N":
                    atom.formal_charge = +1


        if res.resnum == len(seq)+2:
            for atom in res.atom:
                if atom.pdbname.strip() == "N":
                    oldname = atom.pdbname
                    newname = ""
                    for i in oldname:
                        if i == "N":
                            newname += "O"
                        else:
                            newname += i
                    atom.pdbname = newname
                    atom.element = "O"
                    atom.formal_charge = -1
                else:
                    deleteme.append(atom.index)
    st.deleteAtoms(deleteme)
    add_hydrogens(st)

    for chain in st.chain:
        chain.name = "A"
    i = 1
    for res in st.residue:
        res.resnum = i
        i += 1

    return st

##------------------------------------------------------------------------------------------------------------
# For Mod-NZwit peptides: cleans up terminal and creates an acetyl capping group on the c-terminal
def acetyl_peptide(st):
    delete_atoms = []
    add_hydrogen_atoms = []   # update the residue numbers
    end_resnum = len(st.residue)

    # mark ends for protonation/deprotonation
    delete_hydrogens(st)
    add_hydrogens(st)
    for res in st.residue:
        if res.resnum == 1:
            for at in res.atom:
                # print(at.pdbname)
                if at.pdbname == " H1 ":
                    at.element = "C"
                    at.pdbname = " C0 "
        if res.resnum == end_resnum:
            for at in res.atom:
                if at.element == "O" and at.pdbname.strip() == "": # Schrodinger won't print "O1" so I had to use ""
                    at.element = "C"
                    at.pdbname = " C0 "

    delete_hydrogens(st)
    add_hydrogens(st)
    #build.add_hydrogens(st, atom_list=add_hydrogen_atoms)
    #st.deleteAtoms(delete_atoms)
    # st.write("FF_test.pdb")
    return st
    #minimize.minimize_structure(st, cleanup=True)

##------------------------------------------------------------------------------------------------------------
# Clean up the terminal and assign the correct chain name and residue numbers
def zwitter_peptide(st):
    delete_atoms = []
    add_hydrogen_atoms = []   # update the residue numbers
    end_resnum = len(st.residue)

    # mark ends for protonation/deprotonation
    for res in st.residue:
        if res.resnum == 1:
            for at in res.atom:
                if at.pdbname == " N  ":
                    at.formal_charge = 1
                    #add_hydrogen_atoms.append(at.index)
        elif res.resnum == end_resnum:
            for at in res.atom:
                if at.element == "O" and at.pdbname.strip() == "": # Schrodinger won't print "O1" so I had to use ""
                    at.formal_charge = -1

    delete_hydrogens(st)
    add_hydrogens(st)
    #build.add_hydrogens(st, atom_list=add_hydrogen_atoms)
    #st.deleteAtoms(delete_atoms)
    return st
    #minimize.minimize_structure(st, cleanup=True)

##------------------------------------------------------------------------------------------------------------
## peptide structure class
class peptide_struct:
    def __init__(self):
        self.seqs = []
        self.name = ''
        self.Min = False

##------------------------------------------------------------------------------------------------------------
## function to read DNA input file and extract settings/sequence
def file_reader(inputfile):
    inputfile = 'input/' + inputfile
    pep = peptide_struct()

    flag = False

    fh = open(inputfile)
    for line in fh:
        temp = line.split()
        #print(temp)
        if temp[0] == "":
            continue
        if temp[0] == "#":
            if temp[1] == "!Min":
                pep.Min=True
            if temp[1] == "!Name":
                try:
                    pep.name = temp[2]
                except:
                    print("You didn't include a name!")
        if temp[0] == ">":
            flag = True
            continue
        if temp[0] == "<":
            flag = False
            continue
        if flag:
            pep.seqs.append(temp)

    fh.close()
    return pep

##------------------------------------------------------------------------------------------------------------
## function to build beta sheets
def build_sheet(seq):
    st=helix.process_sequence(seq)
    # st = fix_helix(st, seq)
    # st = buildpeptide.build_peptide(seq)
    # st = acetyl_peptide(st)
    # print('hity')
    # st.write("out/VL_preZwitter.mae")
    # st = zwitter_peptide(st)
    return st

#--------------------------------------------
## function to build CARBON CAPPED alpha helix
def build_helix(seq):
    #print(helix.__file__)
    st=helix.process_sequence(seq)
    st = fix_helix(st, seq)

    return st
