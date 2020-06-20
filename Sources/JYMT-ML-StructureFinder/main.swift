//
//  main.swift
//  JYMT-ML-StructureFinder
//
//  Created by Jerry Yan on 3/29/20.
//
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

public let xyz2mol = importXyz2mol()

var (xyzSet, fileName) = xyzFileInput()
print()

var (saveResults, writePath) = exportingPathInput("csv")
print()

let tInitial = Date()
let snt = sntData(from: xyzSet)

if saveResults {
    var csvFile = TextFile()
    csvFile.content = createCSVString(header: ["SMILES", "Validity"], data: snt)
    let csvUrl = writePath.appendingPathComponent("result_" + String(Int(tInitial.timeIntervalSince1970)) + ".csv")
    csvFile.safelyExport(toFile: csvUrl, affix: "csv")
}

let timeTaken = -(Double(tInitial.timeIntervalSinceNow))

print()
print("Computation completed in \(timeTaken.rounded(digitsAfterDecimal: 4)) s.")
print()

print("Exited on \(displayTime(Date())).")
print()
