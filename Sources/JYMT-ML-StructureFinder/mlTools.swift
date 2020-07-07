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
internal typealias CSVData = Dictionary<String, AnyObject>

internal struct SNTTuple {
    internal var smiles: String
    internal var validity: Bool
    
    internal init(_ smiles: String, _ validity: Bool) {
        self.smiles = smiles
        self.validity = validity
    }
}

extension SNTTuple: Hashable {
    internal static func == (lhs: SNTTuple, rhs: SNTTuple) -> Bool {
        return lhs.smiles == rhs.smiles && lhs.validity == rhs.validity
    }
    internal func hash(into hasher: inout Hasher) {
        hasher.combine(smiles)
        hasher.combine(validity)
    }
}

extension SNTTuple {
    internal func convertToCsvData() -> CSVData {
        var csvD = CSVData()
        csvD["SMILES"] = smiles as AnyObject
        csvD["Validity"] = (validity ? 1 : 0) as AnyObject
        return csvD
    }
}

extension SNTTuple {
    internal var csvDataFormat: CSVData {
        return convertToCsvData()
    }
}
