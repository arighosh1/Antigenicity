from quantiprot.metrics.aaindex import get_aa2charge, get_aa2hydropathy
import os
# path = "C:\Users\Arslan\Desktop\binding_affinity\binding_affinity\protein_data.csv"
# os.chdir(path)

import sys
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import re
import streamlit as st
import os



    # basedir = os.path.dirname(os.path.abspath(__file__))
    # basedir = basedir.replace("bioprogram", "media\\file")
    # basedir+='\\'
    # basedir+=sys.argv[1]
    # basedir=re.sub(' /', '//', basedir)
    # aln=open(basedir,'r')
    # aln.read()
file=''
uploaded_files = st.file_uploader("Choose a Protein File To Check Its Antigenicity : ", accept_multiple_files=True)
for uploaded_file in uploaded_files:
    file = uploaded_file.read()

if file!='':

    X = ProteinAnalysis(str(file))
    # print("number of alanine in protein is = ", X.count_amino_acids()['A'])
    # print("number of valine in protein is = ",X.count_amino_acids()['V'])
    # print("number of  isoleucine in protein is = ",X.count_amino_acids()['L'])
    # print("number of leucine in protein is = ",X.count_amino_acids()['I'])
    alanine = X.get_amino_acids_percent()['A']
    valine = X.get_amino_acids_percent()['V']
    isoleucine = X.get_amino_acids_percent()['I']
    leucine = X.get_amino_acids_percent()['L']
    # print(alanine)
    # print(valine)
    # print(isoleucine)
    # print(leucine)
    # print("molecular_weight =", "%0.2f" % X.molecular_weight())
    # print("Aromaticity = ","%0.2f" % X.aromaticity())
    # print("instability_index = ","%0.2f" % X.instability_index())
    # print("Isoelectric point = ","%0.2f" % X.isoelectric_point())
    sec_struc = X.secondary_structure_fraction()
    # print("Secondary structure fraction = ","%0.2f" % sec_struc[0])
    epsilon_prot = X.molar_extinction_coefficient()
    # print(epsilon_prot)

    # Aliphatic index = X(Ala) + a * X(Val) + b * ( X(Ile) + X(Leu) )
    Aliphatic_index = ((alanine + 2.9) * (valine + 3.9) * (isoleucine + leucine)) / 100

    # print(Aliphatic_index)
    index = Aliphatic_index

    # basedir = os.path.dirname(os.path.abspath(__file__))
    # print(basedir)

    if index >= 0.1:
        st.title("The Protein is antigenic protein")
    else:
        st.title("The Protein is non-antigenic protein")

else:
    st.text("Input Your Sequence")