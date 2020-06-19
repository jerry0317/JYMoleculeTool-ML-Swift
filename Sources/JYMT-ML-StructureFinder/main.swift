//
//  main.swift
//  JYMT-ML-StructureFinder
//
//  Created by Jerry Yan on 3/29/20.
//
//  **IMPORTANT**: Set environment variable PYTHON_LIBRARY to the Python path that is able to import RDKit.
//
//  TODO: Automatically search for it?
//


import Foundation
import JYMTBasicKit
import TensorFlow
#if canImport(PythonKit)
    import PythonKit
#else
    import Python
#endif

printWelcomeBanner("Structure Finder w/ ML")

//// Import RDKit
//guard let rdkit = try? Python.attemptImport("rdkit") else {
//    fatalError("RDKit not found. Please set the environment variable PYTHON_LIBRARY to the Python path that is able to import RDKit.")
//}
//
//print("RDKit successfully imported.")

//let X = Tensor<Double>([[1,2,3], [1,5,6], [7,8,9]])
//print(X.count(1))

let sys = Python.import("sys")

guard let xyz2mol = try? Python.attemptImport("xyz2mol") else {
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

var (xyzSet, fileName) = xyzFileInput()

let rawAtoms = xyzSet.atoms!

let combAtoms = rawAtoms.removed(byElement: .hydrogen)

//var (saveResults, writePath) = exportingPathInput("xyz & mol")

let xyzString = XYZFile(fromAtoms: combAtoms).xyzString!
let correctSmiles = String(xyz2mol.xyz2smiles(xyzString))!


//var possibleAtoms = [[rawAtoms[0]]]
//
//for i in 1..<rawAtoms.count {
//    possibleAtoms.append(rawAtoms[i].possibles)
//}

let A1 = selectFarthestAtom(from: combAtoms) ?? rawAtoms[0]
print()

let combrAtoms = combAtoms.removed(A1)
let initialSMol = StrcMolecule(Set([A1]))

var possibleList: [StrcMolecule] = []

let tInitial = Date()
print("Computation started on \(displayTime(tInitial)).")
print()

print("**Structure Filtering Process**")
possibleList = rcsActionDynProgrammed(rAtoms: combrAtoms, stMolList: [initialSMol], filters: [.minimumBondLength, .bondTypeLength, .valence])
//possibleList = rcsActionDynProgrammed(rAtoms: combrAtoms, stMolList: [initialSMol])

let strTimeTaken = -(Double(tInitial.timeIntervalSinceNow))

print()
print("**Structure Filtering completed.**")
print()

var possibleSmiles = [String]()


for psMol in possibleList {
    possibleSmiles.append(String(xyz2mol.xyz2smiles(XYZFile(fromAtoms: Array(psMol.atoms)).xyzString!))!)
}

print()
print("Correct SMILES: \(correctSmiles)")
print("Total number of filtered structures: \(possibleList.count)")
print("Total number of SMILES: \(Set(possibleSmiles).count)")
print("Number of strutures with correct SMILES: \(possibleSmiles.filter({$0 == correctSmiles}).count)")
print("Possible SMILES (first 5): \(Array(Set(possibleSmiles)).prefix(5)) ...")

let smilesTimeTaken = -(Double(tInitial.timeIntervalSinceNow)) - strTimeTaken

print()
print("Structure Filtering completed in \(strTimeTaken.rounded(digitsAfterDecimal: 4)) s.")
print("SMILES computation completed in \(smilesTimeTaken.rounded(digitsAfterDecimal: 4)) s.")
print()

print("Exited on \(displayTime(Date())).")
print()
