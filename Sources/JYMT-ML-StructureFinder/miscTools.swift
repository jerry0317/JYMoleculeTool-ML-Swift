//
//  miscTools.swift
//  JYMT-ML-StructureFinder
//
//  Created by Jerry Yan on 6/20/20.
//

import Foundation
import JYMTBasicKit
#if canImport(PythonKit)
    import PythonKit
#else
    import Python
#endif

func importXyz2mol() -> PythonObject {
    let sys = Python.import("sys")
    guard let xyz2molTemp = try? Python.attemptImport("xyz2mol") else {
        print("Uhh... Failed to import xyz2mol")
        print("Please install xyz2mol in one of the following paths:")
        let paths = sys.path
        for path in paths {
            if String(path)!.contains("site-packages"){
                print(path)
            }
        }
        fatalError()
    }
    return xyz2molTemp
}

func nonHAtoms(from xyzSet: XYZFile) -> [Atom] {
    let rawAtoms = xyzSet.atoms!
    let combAtoms = rawAtoms.removed(byElement: .hydrogen)
    return combAtoms
}

func findCorrectSmiles(from atoms: [Atom]) -> String {
    guard !atoms.isEmpty else {
        return ""
    }
    
    let xyzString = XYZFile(fromAtoms: atoms).xyzString!
    let correctSmiles = String(xyz2mol.xyz2smiles(xyzString))!
    return correctSmiles
}

func findPossibleSmiles(from atoms: [Atom], rcsFilters: Set<StrcFilter> = [.minimumBondLength, .bondTypeLength, .valence]) -> Set<String> {
    guard !atoms.isEmpty else {
        return []
    }
    let A1 = selectFarthestAtom(from: atoms) ?? atoms[0]
    print()

    let combrAtoms = atoms.removed(A1)
    let initialSMol = StrcMolecule(Set([A1]))

    var possibleList: [StrcMolecule] = []
    
    possibleList = rcsActionDynProgrammed(rAtoms: combrAtoms, stMolList: [initialSMol], filters: rcsFilters)
    
    var possibleSmiles = [String]()
    for psMol in possibleList {
        let smiles = String(xyz2mol.xyz2smiles(XYZFile(fromAtoms: Array(psMol.atoms)).xyzString!))!
        possibleSmiles.append(smiles)
    }
    
    return Set(possibleSmiles)
}

func smilesNTruths(uniqueSmiles: Set<String>, correctSmiles: String) -> [Dictionary<String, AnyObject>] {
    var smilesNTruths = [Dictionary<String, AnyObject>]()
    for uSmiles in uniqueSmiles {
        var dict = Dictionary<String, AnyObject>()
        dict["SMILES"] = uSmiles as AnyObject
        dict["Validity"] = (uSmiles == correctSmiles ? 1 : 0) as AnyObject
        smilesNTruths.append(dict)
    }
    return smilesNTruths
}

func sntData(from xyzSet: XYZFile, rcsFilters: Set<StrcFilter> = [.minimumBondLength, .bondTypeLength, .valence]) -> [Dictionary<String, AnyObject>] {
    let atoms = nonHAtoms(from: xyzSet)
    let correctSmiles = findCorrectSmiles(from: atoms)
    let uniqueSmiles = findPossibleSmiles(from: atoms, rcsFilters: rcsFilters)
    return smilesNTruths(uniqueSmiles: uniqueSmiles, correctSmiles: correctSmiles)
}
