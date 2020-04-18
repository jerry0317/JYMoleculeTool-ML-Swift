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

// Import RDKit
guard let rdkit = try? Python.attemptImport("rdkit") else {
    fatalError("RDKit not found. Please set the environment variable PYTHON_LIBRARY to the Python path that is able to import RDKit.")
}

print("RDKit successfully imported.")

let X = Tensor<Double>([[1,2,3], [1,5,6], [7,8,9]])
print(X.count(1))
