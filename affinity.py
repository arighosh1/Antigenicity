from quantiprot.metrics.aaindex import get_aa2charge, get_aa2hydropathy
import os
# path = "C:\Users\Arslan\Desktop\binding_affinity\binding_affinity\protein_data.csv"
# os.chdir(path)
from Bio.SeqUtils.ProtParam import ProteinAnalysis 
X = ProteinAnalysis(input("Enter your sequence: ")) 
# print("number of alanine in protein is = ", X.count_amino_acids()['A'])
# print("number of valine in protein is = ",X.count_amino_acids()['V'])
# print("number of  isoleucine in protein is = ",X.count_amino_acids()['L'])
# print("number of leucine in protein is = ",X.count_amino_acids()['I'])
alanine = X.get_amino_acids_percent()['A']
valine = X.get_amino_acids_percent()['V']
isoleucine = X.get_amino_acids_percent()['I']
leucine = X.get_amino_acids_percent()['L']
print(alanine)
print(valine)
print(isoleucine)
print(leucine)
print("molecular_weight =", "%0.2f" % X.molecular_weight()) 
print("Aromaticity = ","%0.2f" % X.aromaticity()) 
print("instability_index = ","%0.2f" % X.instability_index()) 
print("Isoelectric point = ","%0.2f" % X.isoelectric_point()) 
sec_struc = X.secondary_structure_fraction() 
print("Secondary structure fraction = ","%0.2f" % sec_struc[0]) 
epsilon_prot = X.molar_extinction_coefficient()  
print(epsilon_prot)

#Aliphatic index = X(Ala) + a * X(Val) + b * ( X(Ile) + X(Leu) )
Aliphatic_index = ((alanine + 2.9) * (valine +3.9) * (isoleucine + leucine))/100

print(Aliphatic_index)

if Aliphatic_index >= 0.1:
    print("antigenic protein")
else:
    print("non-antigenic protein")
