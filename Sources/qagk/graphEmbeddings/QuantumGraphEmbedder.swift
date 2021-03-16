import SwiftQuantumComputing
import Foundation
import ComplexModule

#if os(Linux)

import CBLAS_Linux

#else

import Accelerate

#endif

let IDENTITY = try! Matrix(
    [
        [.one, .zero],
        [.zero, .one]
    ]
)

let FLIP = try! Matrix(
    [
        [.zero, .one],
        [.one, .zero]
    ]
)

let P1 = try! Matrix(
    [
        [.zero, .zero],
        [.zero, .one]
    ]
)

let P0 = try! Matrix(
    [
        [.one, .zero],
        [.zero, .zero]
    ]
)

public extension Matrix {
    static func kronekerProduct(lhs: Matrix, rhs: Matrix) -> Matrix { // -> [[Complex]]
        var rows: [[Complex<Double>]] = []
        for l in 0..<lhs.rowCount {
            for k in 0..<rhs.rowCount {
                var row: [Complex<Double>] = []
                for i in 0..<lhs.columnCount {
                    for j in 0..<rhs.columnCount {
                        let lhsCell = lhs[l, i]
                        let rhsCell = rhs[k, j]
                        row.append(Complex<Double>(lhsCell.real * rhsCell.real - lhsCell.imaginary * rhsCell.imaginary, lhsCell.real * rhsCell.imaginary + rhsCell.real * lhsCell.imaginary))
                    }
                }
                rows.append(row)
            }
        }
        return try! Matrix(rows)
    }

    static func multiplyMatrices(lhs: Matrix, rhs: Matrix) -> Matrix { // -> [[Complex]]
        var rows: [[Complex<Double>]] = []
        for l in 0..<lhs.rowCount {
            var row: [Complex<Double>] = []
            for k in 0..<rhs.columnCount {
                let lhsCell = lhs[l, k]
                let rhsCell = rhs[l, k]
                row.append(Complex<Double>(lhsCell.real * rhsCell.real - lhsCell.imaginary * rhsCell.imaginary, lhsCell.real * rhsCell.imaginary + rhsCell.real * lhsCell.imaginary))
            }
            rows.append(row)
        }
        return try! Matrix(rows)
    }

    static func scale(lhs: Matrix, x: Complex<Double>) -> Matrix { // -> [[Complex]]
        var rows: [[Complex<Double>]] = []
        for l in 0..<lhs.rowCount {
            var row: [Complex<Double>] = []
            for k in 0..<lhs.columnCount {
                let lhsCell = lhs[l, k]
                let rhsCell = x
                row.append(Complex<Double>(lhsCell.real * rhsCell.real - lhsCell.imaginary * rhsCell.imaginary, lhsCell.real * rhsCell.imaginary + rhsCell.real * lhsCell.imaginary))
            }
            rows.append(row)
        }
        return try! Matrix(rows)
    }
    
    static func kronekerProduct(matrices: [Matrix]) -> Matrix {
        var product: Matrix = matrices.last!
        for i in (0..<matrices.count - 1).reversed() {
            product = Matrix.kronekerProduct(lhs: matrices[i], rhs: product) 
        }
        return product
    }

    static func multiply(matrices: [Matrix]) -> Matrix {
        var product: Matrix = matrices.first!
        for i in 1..<matrices.count {
            product = Matrix.multiply(lhs: product, rhs: matrices[i])
        }
        return product
    }
}

class ParameterizedGate{
    public var alpha: Double
    public var beta: Double
    public var gamma: Double

    public init() {
        alpha = Double.random(in: 0..<(Double.pi / 2))
        beta = Double.random(in: 0..<(Double.pi / 2))
        gamma = Double.random(in: 0..<(Double.pi / 2))
    }

    public var matrix: Matrix {
        let matrix = try! Matrix(
            [[Complex(cos(beta)*cos(alpha), sin(beta)*cos(alpha)), Complex(cos(gamma)*sin(alpha), sin(gamma)*sin(alpha))],
            [Complex(-cos(gamma)*sin(alpha), sin(gamma)*sin(alpha)), Complex(cos(beta)*cos(alpha), -sin(beta)*cos(alpha))]]
        )
        // print((matrix.krone
        return matrix
    }

    public var asString: String {
        "alpha=\(alpha); beta=\(beta); gamma=\(gamma)"
    }
}

func getControllingQubit(targetQubit: Int, nQubits: Int, offset: Int) -> Int {
    var controlledQubit = targetQubit + offset
    if controlledQubit < 0 {
        controlledQubit = nQubits + controlledQubit
    }
    return controlledQubit
}

func getExtendedMatrix(layers: [[ParameterizedGate]], layer: Int, targetQubit: Int, controllingQubit: Optional<Int> = .none) -> Matrix {
    if let controllingQubitUnwrapped = controllingQubit {
        let p1Part = Matrix.kronekerProduct(
            matrices: (0..<layers[layer].count).map{i -> Matrix in
                if i == targetQubit {
                    return layers[layer][i].matrix
                } else if i == controllingQubitUnwrapped {
                    return P1
                } else {
                    return IDENTITY
                }
            }
        )
        // print(Matrix.kronekerProduct(lhs: P0, rhs: IDENTITY))
        // print(Matrix.kronekerProduct(lhs: IDENTITY, rhs: P0))
        // print(controllingQubitUnwrapped)
        // print("P1")
        // print(p1Part)
        let p0Part = Matrix.kronekerProduct(
            matrices: (0..<layers[layer].count).map{i -> Matrix in
                if i == controllingQubitUnwrapped {
                    return P0
                } else {
                    return IDENTITY
                }
            }
        )
        // print("P0")
        // print(p0Part)
        let result = try! (
            p1Part + p0Part
        ).get()
        // result = Matrix.scale(lhs: result, x: Complex<Double>(1/1.22, 0))
        // print("?????????")
        // print(result)
        return result
    } else {
        return Matrix.kronekerProduct(
            matrices: (0..<layers[layer].count).map{i -> Matrix in
                if i == targetQubit {
                    return layers[layer][i].matrix
                } else {
                    return IDENTITY
                }
            }
        )
    }
}

func makeIdentity(_ nQubits: Int) -> Matrix {
    return try! Matrix(
        (0..<Int(truncating: NSDecimalNumber(decimal: pow(2, nQubits)))).map{i in
            (0..<Int(truncating: NSDecimalNumber(decimal: pow(2, nQubits)))).map{j in
                return j == i ? .one : .zero
            }
        }
    )
}

func getLayerMatrix(layers: [[ParameterizedGate]], layer: Int, offset: Optional<Int> = .none) -> Matrix {
    let nQubits = layers[layer].count
    var product = makeIdentity(nQubits)
    for i in 0..<nQubits {
        var controllingQubit: Optional<Int> = .none
        if let offsetUnwrapped = offset {
            controllingQubit = getControllingQubit(
                targetQubit: i,
                nQubits: nQubits,
                offset: offsetUnwrapped
            )
        }
        let extendedMatrix = getExtendedMatrix(
            layers: layers,
            layer: layer,
            targetQubit: i,
            controllingQubit: controllingQubit
        )
        // print(extendedMatrix)
        product = Matrix.multiply(lhs: product, rhs: extendedMatrix)
        // print("=======================================================================")
        // print(extendedMatrix)
        // print(product)
        // print(Matrix.multiply(lhs: extendedMatrix, rhs: extendedMatrix, rhsTrans: CblasConjTrans))
        // print(Matrix.multiply(lhs: product, rhs: product, rhsTrans: CblasConjTrans))
    }
    // print("=======================================================================")
    // print(Matrix.multiply(lhs: product, rhs: product, rhsTrans: CblasConjTrans))
    return product
}

class QuantumGraphEmbedder {
    // public let graph: Graph
    public let circuit: Circuit
    public let nQubits: Int
    public let gates: [Gate]
    public let predicate: Matrix
    public let parameterizedGates: [[ParameterizedGate]]

    public init(dimensionality nQubits: Int) {
        // Initialize optimized parameters
        // let alpha = Double.random(in: 0..<(Double.pi / 2))
        // let beta = Double.random(in: 0..<(Double.pi / 2))
        // let gamma = Double.random(in: 0..<(Double.pi / 2))
        // let matrix = try! Matrix(
        //     [[Complex(cos(beta)*cos(alpha), sin(beta)*cos(alpha)), Complex(cos(gamma)*sin(alpha), sin(gamma)*sin(alpha))],
        //     [Complex(-cos(gamma)*sin(alpha), sin(gamma)*sin(alpha)), Complex(cos(beta)*cos(alpha), -sin(beta)*cos(alpha))]]
        // )
        var parameterizedGates: [[ParameterizedGate]] = []
        var gates: [Gate] = []
        // First "layer" - apply transforming gates to the qubits representing subjects without conditioning
        func addLayer(offset: Int) {
            var firstLayerGates: [ParameterizedGate] = [] 
            for i in 0...nQubits-1 {
                var gate = ParameterizedGate()
                gates.append(
                    offset == 0 ? .matrix(
                        matrix: gate.matrix, inputs: [i]
                    ) : .controlled(
                        gate: .matrix(
                            matrix: gate.matrix, inputs: [i]
                        ),
                        controls: [
                            getControllingQubit(
                                targetQubit: i,
                                nQubits: nQubits,
                                offset: offset
                            )
                        ]
                    )
                )
                firstLayerGates.append(gate)
                print("Gate applied to the \(i)th qubit: \(gate.asString)")
            }
            parameterizedGates.append(firstLayerGates)
        }
        print("First layer")
        addLayer(offset: 0)
        addLayer(offset: 1)
        addLayer(offset: 2)
        addLayer(offset: 3)
        // print(getExtendedMatrix(layers: parameterizedGates, layer: 0, targetQubit: 1))
        // print(getExtendedMatrix(layers: parameterizedGates, layer: 0, targetQubit: 1, controllingQubit: 0))
        
        // print(getLayerMatrix(layers: parameterizedGates, layer: 0))
        // print(getLayerMatrix(layers: parameterizedGates, layer: 0, offset: -1))
        // print(getLayerMatrix(layers: parameterizedGates, layer: 0, offset: -2))
        // print(getLayerMatrix(layers: parameterizedGates, layer: 0, offset: -3))
        // self.predicate = getLayerMatrix(layers: parameterizedGates, layer: 0)
        self.predicate = Matrix.multiply(
            matrices: [
                getLayerMatrix(layers: parameterizedGates, layer: 0),
                getLayerMatrix(layers: parameterizedGates, layer: 1, offset: -1),
                getLayerMatrix(layers: parameterizedGates, layer: 2, offset: -2),
                getLayerMatrix(layers: parameterizedGates, layer: 3, offset: -3)
            ]
        )
        // var product: Matrix = parameterizedGates.last!.last!.matrix
        // for i in (0..<nQubits - 1).reversed() {

        //     // product = Matrix.kronekerProduct(lhs: parameterizedGates.last![i].matrix, rhs: product) 
        // }
        // print(product)
        // print("Second layer")
        // addLayer(offset: 1)
        // print("Third layer")
        // addLayer(offset: 2)
        // print("Fourth layer")
        // addLayer(offset: 3)
        // print(getControllingQubit(targetQubit: 0, nQubits: 5, offset: 2))
        // self.graph = graph
        // let maxEdgeId = graph.getMaxEdgeId()
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
        self.nQubits = nQubits
        self.parameterizedGates = parameterizedGates
        // print(">>>>>>>>>")
        // print(Matrix.multiply(lhs: predicate, rhs: predicate, rhsTrans: CblasConjTrans))
        // print(Matrix.multiply(lhs: getLayerMatrix(layers: parameterizedGates, layer: 2, offset: -2), rhs: getLayerMatrix(layers: parameterizedGates, layer: 2, offset: -2), rhsTrans: CblasConjTrans))
    }

    func getSubjectMatrix(subject: Int) -> Matrix {
        return makeIdentity(
            self.nQubits
        )
    }

    func getObjectMatrix(object: Int) -> Matrix {
        return makeIdentity(
            self.nQubits
        )
    }

    func run(subject: Int, object: Int) -> CircuitStatevector {
        let U1 = Matrix.multiply(lhs: getSubjectMatrix(subject: subject), rhs: predicate)
        let U2 = getObjectMatrix(object: object)
        let gates: [Gate] = [
            .hadamard(target: 0),
            .controlled(
                gate: .matrix(
                    matrix: U1,
                    inputs: (1...nQubits).map{$0}
                ),
                controls: [0]
            ),
            .not(target: 0),
            .controlled(
                gate: .matrix(
                    matrix: U2,
                    inputs: (1...nQubits).map{$0}
                ),
                controls: [0]
            ),
            .not(target: 0),
            .hadamard(target: 0)
        ]
        let circuit = MainCircuitFactory().makeCircuit(gates: gates)
        return try! circuit.statevector().get()
        // let statevector = try! circuit.statevector().get()
        // return statevector
    }
}
