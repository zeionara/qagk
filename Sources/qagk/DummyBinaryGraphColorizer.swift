import SwiftQuantumComputing

public typealias Graph = Array<(src: Int, dst: Int)>


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


struct DummyBinaryGraphColorizer {
    public let graph: Graph
    public let circuit: Circuit
    public let gates: [Gate]


    public init(_ graph: Graph, inverse: Bool = false) {
        self.graph = graph
        let maxEdgeId = graph.getMaxEdgeId()
        var gates: [Gate] = [
            .hadamard(target: inverse ? maxEdgeId - 0 : 0)
        ]
        for i in 1...graph.getMaxEdgeId() {
            gates.append(.not(target: inverse ? maxEdgeId - i : i))
        }
        for edge in graph {
            gates.append(
                .controlled(
                    gate: .not(target: inverse ? maxEdgeId - edge.dst : edge.dst),
                    controls: [inverse ? maxEdgeId - edge.src : edge.src]
                )
            )
        }
        self.gates = gates
        self.circuit = MainCircuitFactory().makeCircuit(gates: gates)
    }

    public func run() -> CircuitStatevector {
        let statevector = try! circuit.statevector().get()
        // print("Statevector: \(statevector)\n")
        // print("Probabilities: \(statevector.probabilities())\n")
        // print("Summarized probabilities: \(statevector.summarizedProbabilities())\n")
        return statevector
    }
}