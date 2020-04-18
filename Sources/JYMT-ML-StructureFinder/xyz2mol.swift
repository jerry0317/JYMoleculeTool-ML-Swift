//
//  xyz2mol.swift
//  JYMT-ML-StructureFinder
//
//  Created by Jerry Yan on 4/13/20.
//
//  This work is a complete migration from Python to Swift of [xyz2mol](https://github.com/jensengroup/xyz2mol) by Jan H. Jensen.
//
//  The original Implementation was based on the paper:
//
//  Yeonjoon Kim and Woo Youn Kim
//  "Universal Structure Conversion Method for Organic Molecules: From Atomic Connectivity
//  to Three-Dimensional Geometry"
//  Bull. Korean Chem. Soc. 2015, Vol. 36, 1769-1777
//  DOI: 10.1002/bkcs.10334
//
//  Depreacated? (May be invoked later)

import Foundation
import JYMTBasicKit
import TensorFlow

#if canImport(PythonKit)
    import PythonKit
#else
    import Python
#endif

let copy = Python.import("copy")
let itertools = Python.import("itertools")

let rdmolops = Python.import("rdkit.Chem.rdmolops")
let rdEHTTools = try? Python.attemptImport("rdkit.Chem.rdEHTTools")

let defaultdict = Python.import("collections.defaultdict")

let np = Python.import("np")
let networkx = Python.import("nx")

let Chem = Python.import("rdkit.Chem")
let AllChem = Python.import("rdkit.Chem.AllChem")

let __ATOM_LIST__ = ["h",  "he",
"li", "be", "b",  "c",  "n",  "o",  "f",  "ne",
"na", "mg", "al", "si", "p",  "s",  "cl", "ar",
"k",  "ca", "sc", "ti", "v ", "cr", "mn", "fe", "co", "ni", "cu",
"zn", "ga", "ge", "as", "se", "br", "kr",
"rb", "sr", "y",  "zr", "nb", "mo", "tc", "ru", "rh", "pd", "ag",
"cd", "in", "sn", "sb", "te", "i",  "xe",
"cs", "ba", "la", "ce", "pr", "nd", "pm", "sm", "eu", "gd", "tb", "dy",
"ho", "er", "tm", "yb", "lu", "hf", "ta", "w",  "re", "os", "ir", "pt",
"au", "hg", "tl", "pb", "bi", "po", "at", "rn",
"fr", "ra", "ac", "th", "pa", "u",  "np", "pu"]

let atomic_valence: [Int: [Int]] = [
    1: [1],
    6: [4],
    7: [3,4],
    8: [2,1],
    9: [1],
    14: [4],
    15: [5,3], //[5,4,3]
    16: [6,3,2], //[6,4,2]
    17: [1],
    32: [4],
    35: [1],
    53: [1]
]

let atomic_valence_electrons: [Int: Int] = [
    1: 1,
    6: 4,
    7: 5,
    8: 6,
    9: 7,
    14: 4,
    15: 5,
    16: 6,
    17: 7,
    32: 4,
    35: 7,
    53: 7
]

func str_atom(_ atom: Int) -> String {
    return __ATOM_LIST__[atom - 1]
}

func int_atom(_ atom: String) -> Int {
    return __ATOM_LIST__.firstIndex(of: atom.lowercased())! + 1
}

func get_UA(_ maxValence_list: [Int], _ valence_list: [Int]) -> ([Int], [Int]) {
    var UA = [Int]()
    var DU = [Int]()
    for (i, (maxValence, valence)) in zip(maxValence_list, valence_list).enumerated() {
        if (maxValence - valence) <= 0 {
            continue
        }
        UA.append(i)
        DU.append(maxValence - valence)
    }
    return (UA, DU)
}

//func get_BO(_ AC: Matrix, _ UA: [Int], _ DU: [Int], _ valences: [Int], UA_pairs?)

func valences_not_too_large(_ BO: PythonObject, _ valences: [Int]) -> Bool {
    let number_of_bonds_list: [Double] = Array(numpy: BO.sum(axis: 1)) ?? []
    for (valence, number_of_bonds) in zip(valences, number_of_bonds_list){
        if number_of_bonds > Double(valence) {
            return false
        }
    }
    return true
}

// func BO_is_OK

func get_atomic_charge(_ atom: Int, _ atomic_valence_electrons: Int, _ BO_valence: Int) -> Int {
    var charge = 0
    if atom == 1 {
        charge = 1 - BO_valence
    } else if atom == 5 {
        charge = 3 - BO_valence
    } else if atom == 15 && BO_valence == 5 {
        charge = 0
    } else if atom == 16 && BO_valence == 6 {
        charge = 0
    } else {
        charge = atomic_valence_electrons - 8 + BO_valence
    }
    return charge
}

// func clean_charges

//func BO2mol(mol: PythonObject, BO_matrix: PythonObject, atoms: [Int], atomic_valence_electrons)

func set_atomic_charges(_ mol: inout PythonObject, _ atoms: [Int], _ atomic_valence_electrons: [Int: Int], _ BO_valences: [Int], _ BO_Matrix: PythonObject, _ mol_charge: Int) -> PythonObject {
    var q = 0
    for (i, atom) in atoms.enumerated() {
        let a = mol.GetAtomWithIdx(i)
        var charge = get_atomic_charge(atom, atomic_valence_electrons[atom]!, BO_valences[i])
        q += charge
        if atom == 6 {
            let BO_Tensor = Tensor<Double>(numpy: BO_Matrix)!
            let number_of_single_bonds_to_C = BO_Tensor[i, 0...].count(1)
            if number_of_single_bonds_to_C == 2 && BO_valences[i] == 2 {
                q += 1
                charge = 0
            }
            if number_of_single_bonds_to_C == 3 && q + 1 < mol_charge {
                q += 2
                charge = 1
            }
        }
        
        if charge != 0 {
            a.SetFormalCharge(charge)
        }
    }
    
//    mol = clean_charges(mol)
    
    return mol
}

func set_atomic_radicals(_ mol: inout PythonObject, _ atoms: [Int], _ atomic_valence_electrons: [Int: Int], _ BO_valences: [Int]) -> PythonObject {
    for (i, atom) in atoms.enumerated() {
        let a = mol.GetAtomWithIdx(i)
        let charge = get_atomic_charge(atom, atomic_valence_electrons[atom]!, BO_valences[i])
        
        if charge != 0 {
            a.SetNumRadicalElectrons(abs(charge))
        }
    }
    
    return mol
}
