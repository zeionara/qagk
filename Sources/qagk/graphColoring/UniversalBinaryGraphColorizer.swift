import SwiftQuantumComputing


struct UniversalBinaryGraphColorizer: BinaryGraphColorizer {
    public let graph: Graph
    public let circuit: Circuit
    public let gates: [Gate]

    public init(_ graph: Graph) {
        self.graph = graph
        let maxEdgeId = graph.getMaxEdgeId()
        var gates: [Gate] = []
        var hadamarizedQubits = Set<Int>()
        for (i, edge) in graph.enumerated() {
            let edgeQubutId = maxEdgeId + 1 + i
            if !hadamarizedQubits.contains(edge.src) {
                hadamarizedQubits.insert(edge.src)
                gates.append(.hadamard(target: edge.src))
            }
            if !hadamarizedQubits.contains(edge.dst) {
                hadamarizedQubits.insert(edge.dst)
                gates.append(.hadamard(target: edge.dst))
            }
            gates += [
                // If both vertices are of the color with code "1", then set edge qubit to "1"
                .controlled(
                    gate: .not(target: edgeQubutId),
                    controls: [edge.src, edge.dst]
                ),
                // If both vertices are of the color with code "0", then set edge qubit to "1"
                .not(target: edge.src),
                .not(target: edge.dst),
                .controlled(
                    gate: .not(target: edgeQubutId),
                    controls: [edge.src, edge.dst]
                ),
                .not(target: edge.src),
                .not(target: edge.dst),
                // If both vertices are of the same color, then flip one of them
                .controlled(
                    gate: .not(target: edge.dst),
                    controls: [edgeQubutId]
                )
            ]
        }
        self.gates = gates
        self.circuit = MainCircuitFactory().makeCircuit(gates: gates)
    }
}
