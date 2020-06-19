from openeye import oechem

mol = oechem.OEGraphMol()

xstr = """5

C   -2.12929  -2.1542  2.95186
C   -1.84323  -2.94562  4.21691
O   -1.25291  -1.35775  0.84507
O   -3.21948  -2.75475  2.24091
C   -0.9163  -2.09452  2.01966
"""

ims = oechem.oemolistream()
ims.SetFormat(oechem.OEFormat_XYZ)
ims.openstring(xstr)



mols = []
mol = oechem.OEMol()
for mol in ims.GetOEMols():
    mols.append(oechem.OEMol(mol))

oms = oechem.oemolostream()
oms.SetFormat(oechem.OEFormat_SDF)
oms.openstring()

for mol in mols:
    # oechem.OEWriteMolecule(oms, mol)
    print(oechem.OEMolToSmiles(mol))
