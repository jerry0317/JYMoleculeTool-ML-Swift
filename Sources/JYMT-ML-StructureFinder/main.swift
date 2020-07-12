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

printWelcomeBanner("Structure Finder w/ ML [Pre-processing Part]")

internal let xyz2mol = importXyz2mol()

let rcsFilterUsed = rcsFiltersFromCommandLineArg()

var (xyzFiles, fileNames) = xyzFilesInput()
print()

let indexRange = calculationIndexRangeInput(count: xyzFiles.count)
let xyzFilesToUse = Array(xyzFiles[indexRange])
print()

var (_, writePath) = exportingPathInput("csv", isOptional: false)
print()

let numOfThreads = min(sysconf(CInt(_SC_NPROCESSORS_ONLN)) * 2 - 1, xyzFilesToUse.count)

let tInitial = Date()

let sntCsv = sntAction(from: xyzFilesToUse, xyz2mol: xyz2mol, chunkSize: numOfThreads, rcsFilters: rcsFilterUsed)

let csvUrl = writePath.appendingPathComponent("results_\(indexRange == 0...(xyzFiles.count - 1) ? "All" : "\(indexRange.first! + 1)-\(indexRange.last! + 1)")_\(rcsFilterUsed.count)ppf_" + String(Int(tInitial.timeIntervalSince1970)) + ".csv")

exportCsvFile(from: sntCsv, header: ["SMILES", "Validity"], to: csvUrl)

let timeTaken = -(Double(tInitial.timeIntervalSinceNow))

print()
print("Computation completed in \(timeTaken.rounded(digitsAfterDecimal: 4)) s.")
print()

print("Exited on \(displayTime(Date())).")
print()
