// swift-tools-version:5.3
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "qagk",
    dependencies: [
        .package(url: "https://github.com/zeionara/swift-argument-parser.git", .branch("main")),
        .package(url: "https://github.com/zeionara/SwiftQuantumComputing.git", .branch("master")),
        .package(url: "https://github.com/apple/swift-numerics.git", .exact("0.0.8"))
    ],
    targets: [
        // Targets are the basic building blocks of a package. A target can define a module or a test suite.
        // Targets can depend on other targets in this package, and on products in packages this package depends on.
        .target(
            name: "qagk",
            dependencies: [
                .product(name: "ArgumentParser", package: "swift-argument-parser"),
                .product(name: "SwiftQuantumComputing", package: "SwiftQuantumComputing"),
                .product(name: "ComplexModule", package: "swift-numerics")
            ]
        ),
        .testTarget(
            name: "qagkTests",
            dependencies: ["qagk"]),
    ]
)
