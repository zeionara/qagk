import ArgumentParser
import SwiftQuantumComputing
import TensorFlow

#if os(Linux)

import CBLAS_Linux

#else

import Accelerate

#endif

enum AlgorithmError: Error {
    case unsupportedAlgorithm(message: String)
}

struct Test: ParsableCommand {

    @Argument(help: "Name of an algorithm to launch")
    private var algorithm: String = "universal-binary-graph-coloring"

    // @Argument(help: "Path to output file to save generated table")
    // private var outputPath: String

    // @Option(name: .shortAndLong, default: 3, help: "Precision with which values will be printed")
    // var nDecimalPlaces: Int

    private func makeGraphColorizer(_ graph: Graph) throws -> BinaryGraphColorizer {
        if algorithm == "dummy-binary-graph-coloring" {
            return DummyBinaryGraphColorizer(graph)
        } else if algorithm == "universal-binary-graph-coloring" {
            return UniversalBinaryGraphColorizer(graph)
        } else {
            throw AlgorithmError.unsupportedAlgorithm(message: "Algorithm \(algorithm) is unknown to the system")
        }
    }

    mutating func run(_ result: inout [String: Any]) throws {
        print("Running \(self.algorithm) algorithm...")

        // let state = prepareOneQubitState(zeroCoefficient: 0.3, oneCoefficient: 0.7)
        // let state = prepareQubitStates(coefficients: [0.15, 0.13, 0.12, 0.16, 0.11, 0.17, 0.07, 0.09])

        // print(subject.asProbabilities)
        // let state = prepareQubitStates(coefficients: subject.asProbabilities)

        // print(
        //    print(Matrix.multiply(lhs: state, rhs: state, rhsTrans: CblasConjTrans))
        // )

        // let circuit = MainCircuitFactory().makeCircuit(
        //     gates: [
        //         .matrix(
        //             matrix: state, // prepareOneQubitState(zeroCoefficient: 0.5, oneCoefficient: 0.5),
        //             inputs: [0, 1, 2, 3, 4, 5]
        //         )
        //     ]
        // )




        // try! print(circuit.statevector().get().vector)
        // try! print(circuit.statevector().get().summarizedProbabilities())
        let dimensionality = 6
        let embedder = QuantumGraphEmbedder(dimensionality: dimensionality)
        let subject: [Double] = (0..<Int(pow(2.0, Double(dimensionality)))).map{x in Double.random(in: 0..<1)}
        let object: [Double] = (0..<Int(pow(2.0, Double(dimensionality)))).map{x in Double.random(in: 0..<1)}
        let trueLabel = 1.0

        // let gate = ParameterizedGate()
        // let matrix = gate.betaDerivativeMatrix
        // print(Matrix.multiply(lhs: matrix, rhs: matrix, rhsTrans: CblasConjTrans))

        let lr = 0.03
        let nEpochs = 20
        for i in 1..<nEpochs {
            print("Running \(i) epoch...")
            let inferredLabel = embedder.run(
                subject: subject,
                object: object
            ).firstQubitPositiveneStats
            // print(derivatives)
            print("Inferred label \(inferredLabel)")

            for layer in [0, 3] {
                for qubit in 0..<embedder.parameterizedGates[layer].count {
                    let derivatives = embedder.computeDerivatives(
                        subject: subject,
                        object: object,
                        layer: layer,
                        qubit: qubit
                    )
                    embedder.parameterizedGates[layer][qubit].alpha -= (trueLabel - inferredLabel) * lr * derivatives.alpha       
                    embedder.parameterizedGates[layer][qubit].beta -= (trueLabel - inferredLabel) * lr * derivatives.beta
                    embedder.parameterizedGates[layer][qubit].gamma -= (trueLabel - inferredLabel) * lr * derivatives.gamma
                }
            }
        }

        let inferredLabel = embedder.run(
            subject: subject,
            object: object
        ).firstQubitPositiveneStats
        print("Inferred label \(inferredLabel)")

        // .groupedProbabilities(
        //         byQubits: 0..<dimensionality
        //     )
        
        // let identity = try! Matrix(
        //     [[.one, .zero],
        //     [.zero, .one]]
        // )

        // let inverse = try! Matrix(
        //     [[.zero, .one],
        //     [.one, .zero]]
        // )

        // print(Matrix.kronekerProduct(matrices: [P1, inverse]) + Matrix.kronekerProduct(matrices: [P0, identity]))

        // Matrix.kronekerProduct(lhs: identity, rhs: inverse)
        
        // let graph: Graph = [
        //     (0, 1)
        //     // (0, 2),
        //     // (1, 3),
        //     // (2, 4)
        // ]
        // let graph: Graph = [
        //     (0, 1),
        //     (0, 2)
        // ]
        // let colorizer = NGraphColorizer(graph) // try makeGraphColorizer(graph)
        // try print(colorizer.run().groupedProbabilities(byQubits: [0, 1, 2, 3, 4, 5, 14]).get().keys.sorted())
        // print(colorizer.run().summarizedProbabilities())
    }
}

struct Qagk: ParsableCommand {
    static var configuration = CommandConfiguration(
            abstract: "An experimental tool for developing and testing quantum knowledge graph models",
            subcommands: [Test.self],
            defaultSubcommand: Test.self
    )
}

Qagk.main()
