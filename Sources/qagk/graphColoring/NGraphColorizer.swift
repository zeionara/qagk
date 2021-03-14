import SwiftQuantumComputing
import Foundation

public extension ClosedRange where Bound == Int {
    func getAllCombinations() -> Array<(Bound, Bound)> {
        var combinations: Array<(Bound, Bound)> = []
        for i in self {
            if i + 1 <= self.max()! {
                for j in (i + 1)...self.max()! {
                    combinations.append((i, j))
                }
            }
        }
        return combinations
    }
}

struct NGraphColorizer {
    public let graph: Graph
    public let circuit: Circuit
    public let gates: [Gate]

    public init(_ graph: Graph, n: Int = 3) {
        self.graph = graph
        var gates: [Gate] = []
        let maxEdgeId = graph.getMaxEdgeId()
        for vertexIndex in 0...maxEdgeId {
            // let firstIndicatingQubitIndex = (maxEdgeId + 1 + vertexIndex + 1) * n + vertexIndex
            let firstIndicatingQubitIndex = (maxEdgeId + 1 + vertexIndex) * n + vertexIndex
            print("First indicating qubit index is \(firstIndicatingQubitIndex)")
            // Change all qubits to hadamard basis
            for colorIndex in 1...n {
                gates.append(.hadamard(target: vertexIndex * n + colorIndex - 1))
            }
            // Eliminate states in which a vertex is assigned with two colors
            // gates.append(
            //     .controlled(
            //         gate: .not(target: 3),
            //         controls: [0]
            //     )
            // )
            // gates += [
            //     .controlled(
            //         gate: .not(target: 3),
            //         controls: [0]
            //     ),
            //     // .controlled(
            //     //     gate: .not(target: lhsQubitIndex),
            //     //     controls: [firstIndicatingQubitIndex]
            //     // )
            // ]
            // gates += [
            //     .controlled(
            //         gate: .not(target: 3),
            //         controls: [0, 1]
            //     ),
            //     // .controlled(
            //     //     gate: .not(target: lhsQubitIndex),
            //     //     controls: [firstIndicatingQubitIndex]
            //     // )
            // ]
            var nextI = 0
            for (i, (lhs, rhs)) in (1...n).getAllCombinations().enumerated() {
                let lhsQubitIndex = lhs - 1 + vertexIndex * n
                let rhsQubitIndex = rhs - 1 + vertexIndex * n
                gates += [
                    .controlled(
                        gate: .not(target: firstIndicatingQubitIndex + i),
                        controls: [lhsQubitIndex, rhsQubitIndex]
                    ),
                    .controlledNot(
                        target: lhsQubitIndex,
                        control: firstIndicatingQubitIndex + i
                    )
                ]
                nextI = i + 1
            }
            // Eleminate states in which a vertex is assigned with no colors
            let vertexQubitIndices = (1...n).map{$0 - 1 + vertexIndex * n}
            gates += vertexQubitIndices.map{
                .not(target: $0)
            }
            gates += [
                .controlled(
                    gate: .not(target: firstIndicatingQubitIndex + nextI),
                    controls: vertexQubitIndices
                )
            ]
            gates += vertexQubitIndices.map{
                .not(target: $0)
            }
            gates.append(
                .controlled(
                    gate: .not(
                        target: vertexQubitIndices.last!
                    ),
                    controls: [firstIndicatingQubitIndex + nextI]
                )
            )
        }
        // let nextFreeQubut = (maxEdgeId + 1) * 2 * n + maxEdgeId + n
        let nextFreeQubit = (maxEdgeId + 1 + maxEdgeId + 1) * n + maxEdgeId + 1
        print("Next free qubit is \(nextFreeQubit)")
        for (i, edge) in graph.enumerated() {
            for j in 0..<n {
                gates.append(
                    .controlled(
                        gate: .not(
                            target: nextFreeQubit + i
                        ),
                        controls: [edge.src * n + j, edge.dst * n + j]
                    )
                )
            }
            // gates.append(
            //     .controlled(gate: Gate, controls: [Int])
            // )
        }
        // var gates: [Gate] = []
        // var hadamarizedQubits = Set<Int>()
        // for (i, edge) in graph.enumerated() {
        //     let edgeQubutId = maxEdgeId + 1 + i
        //     if !hadamarizedQubits.contains(edge.src) {
        //         hadamarizedQubits.insert(edge.src)
        //         gates.append(.hadamard(target: edge.src))
        //     }
        //     if !hadamarizedQubits.contains(edge.dst) {
        //         hadamarizedQubits.insert(edge.dst)
        //         gates.append(.hadamard(target: edge.dst))
        //     }
        //     gates += [
        //         // If both vertices are of the color with code "1", then set edge qubit to "1"
        //         .controlled(
        //             gate: .not(target: edgeQubutId),
        //             controls: [edge.src, edge.dst]
        //         ),
        //         // If both vertices are of the color with code "0", then set edge qubit to "1"
        //         .not(target: edge.src),
        //         .not(target: edge.dst),
        //         .controlled(
        //             gate: .not(target: edgeQubutId),
        //             controls: [edge.src, edge.dst]
        //         ),
        //         .not(target: edge.src),
        //         .not(target: edge.dst),
        //         // If both vertices are of the same color, then flip one of them
        //         .controlled(
        //             gate: .not(target: edge.dst),
        //             controls: [edgeQubutId]
        //         )
        //     ]
        // }
        self.gates = gates
        self.circuit = MainCircuitFactory().makeCircuit(gates: gates)
    }

    func run() -> CircuitStatevector {
        let statevector = try! circuit.statevector().get()
        return statevector
    }
}
