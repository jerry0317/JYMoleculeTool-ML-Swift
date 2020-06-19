from openeye import oechem

def xyz2smiles(xstr):
    ims = oechem.oemolistream()
    ims.SetFormat(oechem.OEFormat_XYZ)
    ims.openstring(xstr)
    mols = []
    for mol in ims.GetOEMols():
        mols.append(oechem.OEMol(mol))
    if len(mols) == 0:
        return ""
    else:
        mol = oechem.OEMol(mols[0])
        return oechem.OEMolToSmiles(mol)
