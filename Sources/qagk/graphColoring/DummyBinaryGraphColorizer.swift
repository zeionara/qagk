import SwiftQuantumComputing


struct DummyBinaryGraphColorizer: BinaryGraphColorizer {
    public let graph: Graph
    public let circuit: Circuit
    public let gates: [Gate]


    public init(_ graph: Graph) {
        self.graph = graph
        // let maxEdgeId = graph.getMaxEdgeId()
        var gates: [Gate] = [
            .hadamard(target: 0)
        ]
        for i in 1...graph.getMaxEdgeId() {
            gates.append(.not(target: i))
        }
        for edge in graph {
            gates.append(
                .controlled(
                    gate: .not(target: edge.dst),
                    controls: [edge.src]
                )
            )
        }
        self.gates = gates
        self.circuit = MainCircuitFactory().makeCircuit(gates: gates)
    }
}
