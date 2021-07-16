#!/usr/bin/env python
# coding: utf-8

# In[51]:

from quantiprot.metrics.aaindex import get_aa2charge, get_aa2hydropathy
import os
# path = "C:\\Users\\funge\\Documents\\python\\"
# os.chdir(path)
import streamlit as st
from Bio.SeqUtils.ProtParam import ProteinAnalysis
file=''
uploaded_files = st.file_uploader("Choose a Protein File To Check Its Antigenicity : ", accept_multiple_files=True)
for uploaded_file in uploaded_files:
    file = uploaded_file.read()

if file!='':

    X = ProteinAnalysis(str(file))

    # st.write("number of alanine in protein is = ", X.count_amino_acids()['A'])
    # st.write("number of valine in protein is = ",X.count_amino_acids()['V'])
    # st.write("number of  isoleucine in protein is = ",X.count_amino_acids()['L'])
    # st.write("number of leucine in protein is = ",X.count_amino_acids()['I'])
    alanine = X.get_amino_acids_percent()['A']
    valine =  X.get_amino_acids_percent()['V']
    isoleucine =  X.get_amino_acids_percent()['I']
    leucine =  X.get_amino_acids_percent()['L']
    # st.write(alanine)
    # st.write(valine)
    # st.write(isoleucine)
    # st.write(leucine)
    # st.write("molecular_weight =", X.molecular_weight())
    # st.write("Aromaticity = ",X.aromaticity())
    # st.write("instability_index = ",X.instability_index())
    # st.write("Isoelectric point = ",X.isoelectric_point())
    sec_struc = X.secondary_structure_fraction()
    # st.write("Secondary structure fraction = ",sec_struc[0])
    epsilon_prot = X.molar_extinction_coefficient()
    # st.write(epsilon_prot)


    # In[52]:


    #Aliphatic index = X(Ala) + a * X(Val) + b * ( X(Ile) + X(Leu) )
    Aliphatic_index = ((alanine + 2.9) * (valine +3.9) * (isoleucine + leucine))/10


    # In[53]:


    # st.write(Aliphatic_index)


    # In[54]:


    if Aliphatic_index >= 0.1:
        st.title("The Protein is antigenic protein")
    else:
        st.title("The Protein is non-antigenic protein")

else:
    st.text("Input Your Sequence")




