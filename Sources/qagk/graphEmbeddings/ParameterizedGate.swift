import SwiftQuantumComputing
import Foundation
import ComplexModule

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
            [
                [Complex(cos(beta)*cos(alpha), sin(beta)*cos(alpha)), Complex(cos(gamma)*sin(alpha), sin(gamma)*sin(alpha))],
                [Complex(-cos(-gamma)*sin(alpha), -sin(-gamma)*sin(alpha)), Complex(cos(-beta)*cos(alpha), sin(-beta)*cos(alpha))]
            ]
        )
        return matrix
    }

    public var alphaDerivativeMatrix: Matrix {
        let matrix = try! Matrix(
            [
                [Complex(cos(beta)*cos(alpha + Double.pi/2), sin(beta)*cos(alpha + Double.pi/2)), Complex(cos(gamma)*sin(alpha + Double.pi/2), sin(gamma)*sin(alpha + Double.pi/2))],
                [Complex(-cos(-gamma)*sin(alpha + Double.pi/2), -sin(-gamma)*sin(alpha + Double.pi/2)), Complex(cos(-beta)*cos(alpha + Double.pi/2), sin(-beta)*cos(alpha + Double.pi/2))]
            ]
        )
        return matrix
    }

    public func getBetaDerivativeMatrix(_ variant: GateVariant = .normal) -> Matrix {
        func makeSubmatrix(gamma: Double) -> Matrix {
            let matrix = try! Matrix(
                [
                    [Complex(cos(beta + Double.pi/2)*cos(alpha), sin(beta + Double.pi/2)*cos(alpha)), Complex(cos(gamma)*sin(alpha), sin(gamma)*sin(alpha))],
                    [Complex(-cos(-gamma)*sin(alpha), -sin(-gamma)*sin(alpha)), Complex(cos(-(beta + Double.pi/2))*cos(alpha), sin(-(beta + Double.pi/2))*cos(alpha))]
                ]
            )
            return matrix
            // return Matrix.scale(lhs: matrix, x: Complex<Double>(1 / 1.8))
        }
        // return makeSubmatrix(gamma: 0)
        // let unscaledMatrix = try! (makeSubmatrix(gamma: 0) + makeSubmatrix(gamma: Double.pi)).get()
        // let squaredMatrix = Matrix.multiply(lhs: unscaledMatrix, rhs: unscaledMatrix, rhsTrans: CblasConjTrans)
        // print(1 / squaredMatrix[0, 0])
        // let scaledMatrix = Matrix.scale(lhs: unscaledMatrix, x: Complex<Double>(1 / sqrt(squaredMatrix[0, 0].real), squaredMatrix[0, 0].imaginary))
        // return try! (makeSubmatrix(gamma: 0) + makeSubmatrix(gamma: Double.pi)).get().unitary // scaledMatrix
        let result = makeSubmatrix(gamma: variant == .normal ? 0 : Double.pi)
        // print("beta")
        // print(result)
        return result
    }

    public func getGammaDerivativeMatrix(_ variant: GateVariant = .normal) -> Matrix {
        func makeSubmatrix(beta: Double) -> Matrix {
            let matrix = try! Matrix(
                [
                    [Complex(cos(beta)*cos(alpha), sin(beta)*cos(alpha)), Complex(cos(gamma + Double.pi/2)*sin(alpha), sin(gamma + Double.pi/2)*sin(alpha))],
                    [Complex(-cos(-(gamma + Double.pi/2))*sin(alpha), -sin(-(gamma + Double.pi/2))*sin(alpha)), Complex(cos(-beta)*cos(alpha), sin(-beta)*cos(alpha))]
                ]
            )
            return matrix
            // return Matrix.scale(lhs: matrix, x: Complex<Double>(1 / 1.8))
        }
        // return makeSubmatrix(gamma: 0)
        // let unscaledMatrix = try! (makeSubmatrix(gamma: 0) + makeSubmatrix(gamma: Double.pi)).get()
        // let squaredMatrix = Matrix.multiply(lhs: unscaledMatrix, rhs: unscaledMatrix, rhsTrans: CblasConjTrans)
        // print(1 / squaredMatrix[0, 0])
        // let scaledMatrix = Matrix.scale(lhs: unscaledMatrix, x: Complex<Double>(1 / sqrt(squaredMatrix[0, 0].real), squaredMatrix[0, 0].imaginary))
        // return try! (makeSubmatrix(beta: 0) + makeSubmatrix(beta: Double.pi)).get().unitary // scaledMatrix
        let result = makeSubmatrix(beta: variant == .normal ? 0 : Double.pi)
        // print("gamma")
        // print(result)
        return result
    }

    public func getDerivativeMatrix(_ optimizedParameter: OptimizedParameter, variant: GateVariant = .normal) -> Matrix {
        print("Getting derivative matrix...")
        if optimizedParameter == .alpha {
            return alphaDerivativeMatrix
        } else if optimizedParameter == .beta {
            return getBetaDerivativeMatrix(variant)
        } else if optimizedParameter == .gamma {
            return getGammaDerivativeMatrix(variant)
        }
        return matrix
    }

    public var asString: String {
        "alpha=\(alpha); beta=\(beta); gamma=\(gamma)"
    }
}
