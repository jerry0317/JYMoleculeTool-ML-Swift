//
//  mlTools.swift
//  JYMT-ML-StructureFinder
//
//  Created by Jerry Yan on 6/19/20.
//

import Foundation
import JYMTBasicKit
import TensorFlow

func numOfAtoms(_ smiles: String) -> Int {
    return smiles.filter({$0.isUppercase}).count
}

