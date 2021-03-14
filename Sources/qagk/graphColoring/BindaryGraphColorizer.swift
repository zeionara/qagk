import SwiftQuantumComputing

public typealias Graph = Array<(src: Int, dst: Int)>

public protocol BinaryGraphColorizer {
    var graph: Graph {get}
    var circuit: Circuit {get}
    var gates: [Gate] {get}

    init(_ graph: Graph)
}

public extension BinaryGraphColorizer {
    func run() -> CircuitStatevector {
        let statevector = try! circuit.statevector().get()
        return statevector
    }

    func summarize() {
        let statevector = run().groupedProbabilities(byQubits: 0...graph.getMaxEdgeId())
        print(statevector)
    }
}

public extension Graph {
    func getMaxEdgeId() -> Int {
        var result: Optional<Int> = .none
        for edge in self {
            if let unwrappedResult = result {
                if edge.src > unwrappedResult {
                    result = edge.src
                }
                if edge.dst > unwrappedResult {
                    result = edge.dst
                }
            } else {
                if (edge.src >= edge.dst) {
                    result = edge.src
                } else {
                    result = edge.dst
                }
            }
        }
        return result!
    }
}
