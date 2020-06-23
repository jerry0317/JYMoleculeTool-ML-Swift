//
//  mlTools.swift
//  JYMT-ML-StructureFinder
//
//  Created by Jerry Yan on 6/19/20.
//

import Foundation
import JYMTBasicKit
import TensorFlow

internal let supportedElements: Set<ChemElement> = [.hydrogen, .carbon, .oxygen, .nitrogen, .fluorine, .chlorine]
internal typealias SNTData = [Dictionary<String, AnyObject>]
