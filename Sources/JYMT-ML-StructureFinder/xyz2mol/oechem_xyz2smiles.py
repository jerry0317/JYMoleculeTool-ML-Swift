from openeye import oechem, oemolprop

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

    oechem.OEFindRingAtomsAndBonds(mol)
    numInRingAtoms = len(list(mol.GetAtoms(oechem.OEAtomIsInRing())))
    numInRingBonds = len(list(mol.GetBonds(oechem.OEBondIsInRing())))

    # Assign Implicit Hydrogen
    oechem.OEAssignImplicitHydrogens(mol)
    oechem.OEAssignHybridization(mol)
    oechem.OEAssignFormalCharges(mol)

    numOfSingleBonds = 0
    numOfDoubleBonds = 0
    numOfTripleBonds = 0
    numOfQuadrupleBonds = 0

    for bond in mol.GetBonds():
        order = bond.GetOrder()
        if order == 1:
            numOfSingleBonds += 1
        elif order == 2:
            numOfDoubleBonds += 1
        elif order == 3:
            numOfTripleBonds += 1
        elif order == 4:
            numOfQuadrupleBonds += 1

    numsOfAtomWithImplicitHCount = dict()
    numsOfAtomWithDegree = dict()
    numsOfAtomWithExplicitDegree = dict()
    numsOfAtomWithExplicitValence = dict()
    numsOfAtomWithHvyDegree = dict()
    numsOfAtomWithHvyValence = dict()
    numsOfAtomWithValence = dict()

    numsOfAtomWithFormalCharge = dict()

    for k in range(9):
        numsOfAtomWithImplicitHCount[k] = 0
        numsOfAtomWithDegree[k] = 0
        numsOfAtomWithExplicitDegree[k] = 0
        numsOfAtomWithExplicitValence[k] = 0
        numsOfAtomWithHvyDegree[k] = 0
        numsOfAtomWithHvyValence[k] = 0
        numsOfAtomWithValence[k] = 0
        numsOfAtomWithFormalCharge[k - 4] = 0

    numsOfAtomWithHyb = dict()

    for k in range(6):
        numsOfAtomWithHyb[k] = 0

    numOfCAtoms = 0
    numOfNonCHAtoms = 0

    for atom in mol.GetAtoms():
        numsOfAtomWithImplicitHCount[atom.GetImplicitHCount()] += 1
        numsOfAtomWithDegree[atom.GetDegree()] += 1
        numsOfAtomWithExplicitDegree[atom.GetExplicitDegree()] += 1
        numsOfAtomWithExplicitValence[atom.GetExplicitValence()] += 1
        numsOfAtomWithHvyDegree[atom.GetHvyDegree()] += 1
        numsOfAtomWithHvyValence[atom.GetHvyValence()] += 1
        numsOfAtomWithValence[atom.GetValence()] += 1
        numsOfAtomWithHyb[atom.GetHyb()] += 1
        numsOfAtomWithFormalCharge[atom.GetFormalCharge()] += 1
        atomicNum = atom.GetAtomicNum()
        if atomicNum == 6:
            numOfCAtoms += 1
        elif atomicNum != 1:
            numOfNonCHAtoms += 1

    mol2dPSA = oemolprop.OEGet2dPSA(mol)
    molAnionicCarbonCount = oemolprop.OEGetAnionicCarbonCount(mol)
    molAromaticRingCount = oemolprop.OEGetAromaticRingCount(mol)
    # molFractionCsp3 = oemolprop.OEGetFractionCsp3(mol)
    molHalideFraction = oemolprop.OEGetHalideFraction(mol)
    molHBondAcceptorCount = oemolprop.OEGetHBondAcceptorCount(mol)
    molHBondDonorCount = oemolprop.OEGetHBondDonorCount(mol)
    molLipinskiAcceptorCount = oemolprop.OEGetLipinskiAcceptorCount(mol)
    molLipinskiDonorCount = oemolprop.OEGetLipinskiDonorCount(mol)
    molLongestUnbranchedHeavyAtomsChain = oemolprop.OEGetLongestUnbranchedHeavyAtomsChain(mol)
    molLongestUnbranchedCarbonsChain = oemolprop.OEGetLongestUnbranchedCarbonsChain(mol)
    molNumUnspecifiedAtomStereos = oemolprop.OEGetNumUnspecifiedAtomStereos(mol)
    molNumUnspecifiedBondStereos = oemolprop.OEGetNumUnspecifiedBondStereos(mol)

    numOfRotorBonds = oemolprop.OEGetRotatableBondCount(mol)

    result = {
        "SMILES": smiles,
        # "dimension": dimension,
        "numAtoms": numAtoms,
        "numBonds": numBonds,
        "numGroups": numGroups,
        "numAromaticAtoms": numAromaticAtoms,
        "numAromaticBonds": numAromaticBonds,
        "numInRingAtoms": numInRingAtoms,
        "numInRingBonds": numInRingBonds,
        "numOfSingleBonds": numOfSingleBonds,
        "numOfDoubleBonds": numOfDoubleBonds,
        "numOfTripleBonds": numOfTripleBonds,
        "numOfQuadrupleBonds": numOfQuadrupleBonds,
        "numOfCAtoms": numOfCAtoms,
        "numOfNonCHAtoms": numOfNonCHAtoms,
        "numOfRotorBonds": numOfRotorBonds,
        "mol2dPSA": mol2dPSA,
        "molAnionicCarbonCount": molAnionicCarbonCount,
        "molAromaticRingCount": molAromaticRingCount,
        # "molFractionCsp3": molFractionCsp3,
        "molHalideFraction": molHalideFraction,
        "molHBondAcceptorCount": molHBondAcceptorCount,
        "molHBondDonorCount": molHBondDonorCount,
        "molLipinskiAcceptorCount": molLipinskiAcceptorCount,
        "molLipinskiDonorCount": molLipinskiDonorCount,
        "molLongestUnbranchedHeavyAtomsChain": molLongestUnbranchedHeavyAtomsChain,
        "molLongestUnbranchedCarbonsChain": molLongestUnbranchedCarbonsChain,
        "molNumUnspecifiedAtomStereos": molNumUnspecifiedAtomStereos,
        "molNumUnspecifiedBondStereos": molNumUnspecifiedBondStereos
    }

    for k in range(9):
        result["numsOfAtomWithImplicitHCount" + str(k)] = numsOfAtomWithImplicitHCount[k]
        result["numsOfAtomWithDegree" + str(k)] = numsOfAtomWithDegree[k]
        result["numsOfAtomWithExplicitDegree" + str(k)] = numsOfAtomWithExplicitDegree[k]
        result["numsOfAtomWithExplicitValence" + str(k)] = numsOfAtomWithExplicitValence[k]
        result["numsOfAtomWithHvyDegree" + str(k)] = numsOfAtomWithHvyDegree[k]
        result["numsOfAtomWithHvyValence" + str(k)] = numsOfAtomWithHvyValence[k]
        result["numsOfAtomWithValence" + str(k)] = numsOfAtomWithValence[k]
        result["numsOfAtomWithFormalCharge" + str(k - 4)] = numsOfAtomWithFormalCharge[k - 4]

    for k in range(6):
        result["numsOfAtomWithHyb" + str(k)] = numsOfAtomWithHyb[k]

    return result

def smiles2features(sstr):
    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, sstr)
    return molobj2features(mol)

def nullifyOEThrowStream():
    oechem.OEThrow.SetOutputStream(oechem.oenul)
