// swift-tools-version:5.1
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "JYMoleculeTool-ML-Swift",
    platforms: [
        .macOS(.v10_15)
    ],
    products: [
        .executable(name: "JYMT-ML-StructureFinder", targets: ["JYMT-ML-StructureFinder"])
    ],
    dependencies: [
        .package(url: "https://github.com/jerry0317/JYMTKit", .branch("master"))
    ],
    targets: [
        .target(
            name: "JYMT-ML-StructureFinder",
            dependencies: ["JYMTBasicKit"])
    ]
)
