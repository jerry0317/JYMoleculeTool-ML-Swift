//
//  tensorFlowExtensions.swift
//  JYMT-ML-StructureFinder
//
//  Created by Jerry Yan on 4/18/20.
//

import Foundation
import TensorFlow

extension Tensor where Scalar : Equatable {
    func count(_ n: Scalar) -> Int{
        return(self.gathering(where: self .== n)).scalarCount
    }
}
