from openeye import oechem

# Find the molecule object based on the .xyz input. Return None if nothing found.
def xyz2molobj(xstr):
    ims = oechem.oemolistream()
    ims.SetFormat(oechem.OEFormat_XYZ)
    ims.openstring(xstr)
    mols = []
    for mol in ims.GetOEMols():
        mols.append(oechem.OEMol(mol))
    if len(mols) == 0:
        return None
    else:
        return oechem.OEMol(mols[0])

# Find the SMILES string based on the .xyz input. Return empty string if nothing found.
def xyz2smiles(xstr):
    mol = xyz2molobj(xstr)
    if mol == None:
        return ""
    else:
        return oechem.OEMolToSmiles(mol)

# Find related features from molobj.
def molobj2features(mol):
    smiles = oechem.OEMolToSmiles(mol)
    # dimension = oechem.OEGetDimensionFromCoords(mol)
    numAtoms = mol.NumAtoms()
    numBonds = mol.NumBonds()
    numGroups = mol.NumGroups()
    oechem.OEAssignAromaticFlags(mol)
    numAromaticAtoms = len(list(mol.GetAtoms(oechem.OEIsAromaticAtom())))
    numAromaticBonds = len(list(mol.GetBonds(oechem.OEIsAromaticBond())))

    result = {
        "SMILES": smiles,
        # "dimension": dimension,
        "numAtoms": numAtoms,
        "numBonds": numBonds,
        "numGroups": numGroups,
        "numAromaticAtoms": numAromaticAtoms,
        "numAromaticBonds": numAromaticBonds
    }

    return result

def smiles2features(sstr):
    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, sstr)
    return molobj2features(mol)
