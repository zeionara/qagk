import ArgumentParser
import SwiftQuantumComputing
import TensorFlow

struct Test: ParsableCommand {

    @Argument(help: "Name of an algorithm to launch")
    private var algorithm: String = "dummy-binary-graph-coloring"

    // @Argument(help: "Path to output file to save generated table")
    // private var outputPath: String

    // @Option(name: .shortAndLong, default: 3, help: "Precision with which values will be printed")
    // var nDecimalPlaces: Int

    mutating func run(_ result: inout [String: Any]) throws {
        print("Running \(self.algorithm) algorithm...")
        let graph: Graph = [
            (0, 1),
            (0, 2),
            (1, 3),
            (2, 4)
        ]
        let colorizer = DummyBinaryGraphColorizer(graph)
        print(colorizer.run().summarizedProbabilities())
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
