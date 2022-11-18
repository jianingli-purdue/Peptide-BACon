import dev
import box
import sys

# script to create pdb and mae files of peptides

##------------------------------------------------------------------------------------------------------------
# peptide variables
seq = "FF"
shape = "beta"
name = "FF_single"

# generate and save
if shape == "alpha":
    st = dev.build_helix(seq)
    st.write(name + ".mae")
    st.write(name + ".pdb")
elif shape == "beta":
    st = dev.build_sheet(seq)
    st.write(name + ".mae")
    st.write(name + ".pdb")
