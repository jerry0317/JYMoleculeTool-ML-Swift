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

func smilesNTruths(uniqueSmiles: Set<String>, correctSmiles: String) -> Set<SNTTuple> {
    var smilesNTruths = Set<SNTTuple>()
    for uSmiles in uniqueSmiles {
        let snt = SNTTuple(uSmiles, uSmiles == correctSmiles)
        smilesNTruths.insert(snt)
    }
    return smilesNTruths
}

func sntAction(from xyzSets: [XYZFile], xyz2mol: PythonObject, chunkSize: Int = 4, rcsFilters: Set<StrcFilter> = [.minimumBondLength, .bondTypeLength, .valence]) -> [CSVData] {
    var snt = Set<SNTTuple>()
    var correctSmilesSet = Set<String>()
    let serialQueue = DispatchQueue(label: "SerialQueue")
    let xyzChunked = xyzSets.shuffled().chunked(by: chunkSize)
    var numOfInValidXyz = 0
    var numOfValidXyz = 0
    var numOfEmptyStrcs = 0
    
    let tInitial = Date()
    
    func printSntProgress() {
        #if DEBUG
        #else
        let numOfXyzProcessed = numOfValidXyz + numOfInValidXyz + numOfEmptyStrcs
        let portionCompleted = Double(numOfXyzProcessed) / Double(xyzSets.count)
        let timeElapsed = -Double(tInitial.timeIntervalSinceNow)
        let eta = ((1 - portionCompleted) / portionCompleted) * timeElapsed
        printStringInLine(toPrintWithSpace("Processing: \(numOfXyzProcessed)/\(xyzSets.count) - \((portionCompleted * 100).rounded(digitsAfterDecimal: 1)) % - ETA: \(timeIntervalToString(eta, maximumUnitCount: 2))", 64))
        #endif
    }
    
    xyz2mol.nullifyOEThrowStream()
    
    DispatchQueue.concurrentPerform(iterations: xyzChunked.count, execute: { i in
        var xyzChunk = [XYZFile]()
        serialQueue.sync {
            xyzChunk = xyzChunked[i]
        }
        
        for xyzSet in xyzChunk {
            let atoms = nonHAtoms(from: xyzSet)
            guard Set(atoms.map({$0.element})).isSubset(of: supportedElements) else {
                serialQueue.sync {
                    numOfInValidXyz += 1
                    printSntProgress()
                }
//                print("Found unsupported element. Data set neglected.")
                continue
            }
            var correctSmiles = ""
            var uniqueSmiles = Set<String>()
            var gCache = GlobalCache()
            let possibleStructures = findPossibleStructures(from: atoms, cache: &gCache, rcsFilters: rcsFilters)
            serialQueue.sync {
                if possibleStructures.isEmpty {
                    numOfEmptyStrcs += 1
                } else {
                    correctSmiles = findCorrectSmiles(from: atoms, xyz2mol: xyz2mol)
                    correctSmilesSet.insert(correctSmiles)
                    uniqueSmiles = findPossibleSmiles(from: possibleStructures, xyz2mol: xyz2mol)
                    snt.formUnion(smilesNTruths(uniqueSmiles: uniqueSmiles, correctSmiles: correctSmiles))
                    numOfValidXyz += 1
                }
                printSntProgress()
            }
        }
    })
    
    print("\n")
    print("Number of input XYZ files: \(numOfValidXyz + numOfInValidXyz + numOfEmptyStrcs)")
    print("Number of valid XYZ files: \(numOfValidXyz)")
    print("Number of unique compounds: \(correctSmilesSet.count)")
    print("Number of possible SMILES: \(snt.count)")
    print()
    
    return snt.map({$0.csvDataFormat})
}

func exportCsvFile(from snt: [CSVData], to csvUrl: URL) {
    var csvFile = TextFile()
    csvFile.content = createCSVString(header: ["SMILES", "Validity"], data: snt)
    csvFile.safelyExport(toFile: csvUrl, affix: "csv")
}

func calculationIndexRangeInput(count: Int) -> ClosedRange<Int> {
    var pass = false
    var result = 0...(count - 1)
    while !pass {
        let inStr = input(name: "calculation range", type: "string", defaultValue: "1-\(count)", printAfterSec: false)
        if inStr.isEmpty {
            pass = true
        } else if !inStr.contains("-") {
            pass = false
            print("Illegal format. Please try again.")
        } else {
            let strSplit = inStr.split(separator: "-", maxSplits: 1, omittingEmptySubsequences: false)
            let pstr = strSplit.map({$0.trimmingCharacters(in: .whitespaces)})
            var startIndex: Int? = nil
            var endIndex: Int? = nil
            if pstr[0].isEmpty {
                startIndex = 0
            } else {
                startIndex = Int(pstr[0])
            }
            if pstr[1].isEmpty {
                endIndex = count - 1
            } else {
                endIndex = Int(pstr[1])
            }
            
            if (startIndex == nil || endIndex == nil) {
                pass = false
                print("Illegal format. Please try again.")
            } else if (startIndex! <= 0 || endIndex! > count || endIndex! < startIndex!) {
                pass = false
                print("Incorrect range. Please try again.")
            } else {
                result = (startIndex! - 1)...(endIndex! - 1)
                pass = true
            }
        }
    }
    print("The calculation range is set to \(result.first! + 1)-\(result.last! + 1).")
    return result
}
