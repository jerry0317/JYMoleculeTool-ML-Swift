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

func findCorrectSmiles(from atoms: [Atom], xyz2mol: PythonObject) -> String {
    guard !atoms.isEmpty else {
        return ""
    }
    
    let xyzString = XYZFile(fromAtoms: atoms).xyzString!
    let correctSmiles = String(xyz2mol.xyz2smiles(xyzString))!
    return correctSmiles
}

func findPossibleStructures(from atoms: [Atom], cache: inout GlobalCache, rcsFilters: Set<StrcFilter> = [.minimumBondLength, .bondTypeLength, .valence], toPrint: Bool = false) -> [StrcMolecule] {
    guard !atoms.isEmpty else {
        return []
    }
    let A1 = selectFarthestAtom(from: atoms) ?? atoms[0]
    if toPrint {
        print()
    }

    let combrAtoms = atoms.removed(A1)
    let initialSMol = StrcMolecule(Set([A1]))

    var possibleList: [StrcMolecule] = []
    
    possibleList = rcsActionDynProgrammed(rAtoms: combrAtoms, stMolList: [initialSMol], filters: rcsFilters, toPrint: toPrint, cache: &cache)
    
    return possibleList
}

func findPossibleSmiles(from strcMolecules: [StrcMolecule], xyz2mol: PythonObject) -> Set<String> {
    var possibleSmiles = [String]()
    for psMol in strcMolecules {
        let smiles = String(xyz2mol.xyz2smiles(XYZFile(fromAtoms: Array(psMol.atoms)).xyzString!))!
        possibleSmiles.append(smiles)
    }
    
    return Set(possibleSmiles)
}

func smilesNTruths(uniqueSmiles: Set<String>, correctSmiles: String) -> [SNTData] {
    var smilesNTruths = [SNTData]()
    for uSmiles in uniqueSmiles {
        var dict = SNTData()
        dict["SMILES"] = uSmiles as AnyObject
        dict["Validity"] = (uSmiles == correctSmiles ? 1 : 0) as AnyObject
        smilesNTruths.append(dict)
    }
    return smilesNTruths
}

func sntAction(from xyzSets: [XYZFile], xyz2mol: PythonObject, chunkSize: Int = 4, rcsFilters: Set<StrcFilter> = [.minimumBondLength, .bondTypeLength, .valence]) -> [SNTData] {
    var snt = [SNTData]()
    let serialQueue = DispatchQueue(label: "SerialQueue")
    let xyzChunked = xyzSets.chunked(by: chunkSize)
    var numOfXyzProcessed = 0
    var numOfValidXyz = 0
    
    let tInitial = Date()
    
    xyz2mol.nullifyOEThrowStream()
    
    DispatchQueue.concurrentPerform(iterations: xyzChunked.count, execute: { i in
        var xyzChunk = [XYZFile]()
        serialQueue.sync {
            xyzChunk = xyzChunked[i]
        }
        
        for xyzSet in xyzChunk {
            let atoms = nonHAtoms(from: xyzSet)
            serialQueue.sync {
                numOfXyzProcessed += 1
            }
            guard Set(atoms.map({$0.element})).isSubset(of: supportedElements) else {
//                print("Found unsupported element. Data set neglected.")
                continue
            }
            var correctSmiles = ""
            var uniqueSmiles = Set<String>()
            var gCache = GlobalCache()
            let possibleStructures = findPossibleStructures(from: atoms, cache: &gCache, rcsFilters: rcsFilters)
            serialQueue.sync {
                correctSmiles = findCorrectSmiles(from: atoms, xyz2mol: xyz2mol)
                uniqueSmiles = findPossibleSmiles(from: possibleStructures, xyz2mol: xyz2mol)
                snt.append(contentsOf: smilesNTruths(uniqueSmiles: uniqueSmiles, correctSmiles: correctSmiles))
                numOfValidXyz += 1
                
                #if DEBUG
                #else
                let portionCompleted = Double(numOfXyzProcessed) / Double(xyzSets.count)
                let timeElapsed = -Double(tInitial.timeIntervalSinceNow)
                let eta = ((1 - portionCompleted) / portionCompleted) * timeElapsed
                printStringInLine(toPrintWithSpace("Calculating: \(numOfXyzProcessed)/\(xyzSets.count) - \((portionCompleted * 100).rounded(digitsAfterDecimal: 1)) % - ETA: \(Int(eta)) s", 70))
                #endif
            }
        }
    })
    
    printStringInLine("")
    print("Number of valid XYZ files: \(numOfValidXyz)")
    print("Number of SMILES calculated: \(snt.count)")
    print()
    
    return snt
}

func exportCsvFile(from snt: [SNTData], to csvUrl: URL) {
    var csvFile = TextFile()
    csvFile.content = createCSVString(header: ["SMILES", "Validity"], data: snt)
    csvFile.safelyExport(toFile: csvUrl, affix: "csv")
}
