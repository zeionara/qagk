import SwiftQuantumComputing
import ComplexModule
import Foundation

func prepareOneQubitState(zeroCoefficient: Double, oneCoefficient: Double) -> Matrix {
    let n = sqrt(zeroCoefficient / oneCoefficient)
    let z = sqrt(1.0 / (n*n + 1))
    let p = sqrt(1.0 - z*z)
    let y = sqrt(1.0 - n*n*z*z)
    let x = sqrt(1.0 - y*y)

    // print("x=\(x), y=\(y), z=\(z), p=\(p), n=\(n)")

    return try! Matrix(
        [
            [Complex<Double>(-x, 0), Complex<Double>(y, 0)],
            [Complex<Double>(-z, 0), Complex<Double>(-p, 0)]
        ]
    )
}

func prepareQubitStates(coefficients: [Double]) -> Matrix {
    if (coefficients.count == 2) {
        return prepareOneQubitState(zeroCoefficient: coefficients[0], oneCoefficient: coefficients[1])
    } else {
        let boundary = coefficients.count / 2
        let zeroStateMatrix = prepareQubitStates(coefficients: Array(coefficients[0..<boundary]))
        let oneStateMatrix = prepareQubitStates(coefficients: Array(coefficients[boundary..<coefficients.count]))
        let firstQubit = prepareOneQubitState(
            zeroCoefficient: Array(coefficients[0..<boundary]).reduce(0, {x, y in x + pow(y, 1)}),
            oneCoefficient: Array(coefficients[boundary..<coefficients.count]).reduce(0, {x, y in x + pow(y, 1)})
        )
        let identity = makeIdentity(Int(log2(Double(boundary))))

        // Construct a state preparation matrix for cases in which current qubit is one
        let oneStateP1 = Matrix.kronekerProduct(
            matrices: [
                P1,
                oneStateMatrix
            ]
        )
        
        let oneStateP0 = Matrix.kronekerProduct(
            matrices: [
                P0,
                identity
            ]
        )
        let oneStateConditionedMatrix = try! (
            oneStateP1 + oneStateP0
        ).get()
        
        // Construct a state preparation matrix for cases in which current qubit is zero
        let inverse = Matrix.kronekerProduct(lhs: FLIP, rhs: identity)
        let zeroStateP1 = Matrix.kronekerProduct(
            matrices: [
                P1,
                zeroStateMatrix
            ]
        )
        let zeroStateP0 = Matrix.kronekerProduct(
            matrices: [
                P0,
                identity
            ]
        )
        let zeroStateConditionedMatrix = Matrix.multiply(
            matrices: [
                inverse,
                try! (
                    zeroStateP1 + zeroStateP0
                ).get(),
                inverse
            ]
        )
        return Matrix.stackAsGates(
            [
                Matrix.kronekerProduct(matrices: [firstQubit, identity]),
                oneStateConditionedMatrix,
                zeroStateConditionedMatrix
            ]
        )
    }
}
