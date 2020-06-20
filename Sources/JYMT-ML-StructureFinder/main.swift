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
print()

var (saveResults, writePath) = exportingPathInput("csv")
print()

let rawAtoms = xyzSet.atoms!

let combAtoms = rawAtoms.removed(byElement: .hydrogen)

let xyzString = XYZFile(fromAtoms: combAtoms).xyzString!
let correctSmiles = String(xyz2mol.xyz2smiles(xyzString))!

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
var smilesNTruths = [Dictionary<String, AnyObject>]()


for psMol in possibleList {
    let smiles = String(xyz2mol.xyz2smiles(XYZFile(fromAtoms: Array(psMol.atoms)).xyzString!))!
    possibleSmiles.append(smiles)
}

let uniqueSmiles = Set(possibleSmiles)

for uSmiles in uniqueSmiles {
    var dict = Dictionary<String, AnyObject>()
    dict["SMILES"] = uSmiles as AnyObject
    dict["Validity"] = (uSmiles == correctSmiles ? 1 : 0) as AnyObject
    smilesNTruths.append(dict)
}

let header = ["SMILES", "Validity"]

func createCSVString(header: [String], data: [Dictionary<String, AnyObject>], nilString: String = "N/A") -> String {
    var csvStr = header.joined(separator: ",") + "\n"
    for dict in data {
        csvStr += header.map({String(describing: dict[$0] ?? nilString as AnyObject)}).joined(separator: ",") + "\n"
    }
    return csvStr
}

if saveResults {
    var csvFile = TextFile()
    csvFile.content = createCSVString(header: header, data: smilesNTruths)
    let csvUrl = writePath.appendingPathComponent("result_" + String(Int(tInitial.timeIntervalSince1970)) + ".csv")
    csvFile.safelyExport(toFile: csvUrl)
}

//print(createCSVString(header: header, data: smilesNTruths))

print()
print("Correct SMILES: \(correctSmiles)")
print("Total number of filtered structures: \(possibleList.count)")
print("Total number of SMILES: \(uniqueSmiles.count)")
print("Number of strutures with correct SMILES: \(possibleSmiles.filter({$0 == correctSmiles}).count)")
print("Possible SMILES (first 5): \(Array(uniqueSmiles).prefix(5)) ...")

let smilesTimeTaken = -(Double(tInitial.timeIntervalSinceNow)) - strTimeTaken

print()
print("Structure Filtering completed in \(strTimeTaken.rounded(digitsAfterDecimal: 4)) s.")
print("SMILES computation completed in \(smilesTimeTaken.rounded(digitsAfterDecimal: 4)) s.")
print()

print("Exited on \(displayTime(Date())).")
print()
