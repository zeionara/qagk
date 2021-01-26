import ArgumentParser
import SwiftQuantumComputing

struct Test: ParsableCommand {

    // @Argument(help: "Path to file containing input report")
    // private var inputPath: String

    // @Argument(help: "Path to output file to save generated table")
    // private var outputPath: String

    // @Option(name: .shortAndLong, default: 3, help: "Precision with which values will be printed")
    // var nDecimalPlaces: Int

    mutating func run(_ result: inout [String: Any]) throws {
        let matrix = try! Matrix(
            [[.one, .zero, .zero, .zero],
            [.zero, .one, .zero, .zero],	
            [.zero, .zero, .zero, .one],
            [.zero, .zero, .one, .zero]]
        )

        let gates: [Gate] = [
            .not(target: 0),
            .hadamard(target: 1),
            .phaseShift(radians: 0.25, target: 2),
            .rotation(axis: .z, radians: 1, target: 3),
            .matrix(matrix: matrix, inputs: [3, 2]),
            .matrix(matrix: matrix, inputs: [0, 3]),
            .oracle(truthTable: ["01", "10"],
                    controls: [0, 1],
                    gate: .rotation(axis: .x, radians: 0.5, target: 3)),
            .oracle(truthTable: ["0"], controls: [0], target: 2),
            .controlled(gate: .hadamard(target: 4), controls: [2]),
            .controlled(gate: .matrix(matrix: matrix, inputs: [4, 2]), controls: [1, 0]),
            .controlledNot(target: 0, control: 3)
        ]

        let circuit = MainCircuitFactory().makeCircuit(gates: gates)

        let statevector = try! circuit.statevector().get()
        print("Statevector: \(statevector)\n")
        print("Probabilities: \(statevector.probabilities())\n")
        print("Summarized probabilities: \(statevector.summarizedProbabilities())\n")
        let groupedProbs = try! statevector.groupedProbabilities(byQubits: [1, 0],
                                                            summarizedByQubits: [4, 3, 2]).get()
        print("Grouped probabilities: \(groupedProbs)")
        print("Unitary: \(try! circuit.unitary().get())\n")
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
