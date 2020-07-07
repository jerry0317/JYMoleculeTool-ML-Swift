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

internal let xyz2mol = importXyz2mol()

var rcsFilterUsed: Set<StrcFilter> = [.minimumBondLength, .bondTypeLength, .valence, .bondAngle, .coplanarity]

if CommandLine.arguments.count >= 2 {
    if CommandLine.arguments.contains("-3ppf") {
        rcsFilterUsed = [.minimumBondLength, .bondTypeLength, .valence]
        print("[Running with 3 pre-processed filters]\n")
    } else if CommandLine.arguments.contains("-5ppf") {
        rcsFilterUsed = [.minimumBondLength, .bondTypeLength, .valence, .bondAngle, .coplanarity]
        print("[Running with 5 pre-processed filters]\n")
    }
}

let numOfThreads = sysconf(CInt(_SC_NPROCESSORS_ONLN)) * 2 - 1

var (xyzFiles, fileNames) = xyzFilesInput()
print()

var (_, writePath) = exportingPathInput("csv", isOptional: false)
print()

print("Note: ~\(numOfThreads) threads will be used for calculation.")
print()

let tInitial = Date()

let sntCsv = sntAction(from: xyzFiles, xyz2mol: xyz2mol, chunkSize: numOfThreads, rcsFilters: rcsFilterUsed)
let csvUrl = writePath.appendingPathComponent("results_\(rcsFilterUsed.count)ppf_" + String(Int(tInitial.timeIntervalSince1970)) + ".csv")

exportCsvFile(from: sntCsv, to: csvUrl)

let timeTaken = -(Double(tInitial.timeIntervalSinceNow))

print()
print("Computation completed in \(timeTaken.rounded(digitsAfterDecimal: 4)) s.")
print()

print("Exited on \(displayTime(Date())).")
print()
