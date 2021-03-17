import SwiftQuantumComputing
import ComplexModule
import Foundation

func prepareOneQubitState(zeroCoefficient: Double, oneCoefficient: Double) -> Matrix {
    print(zeroCoefficient, oneCoefficient)

    let n = sqrt(zeroCoefficient / oneCoefficient)
    let z = sqrt(1.0 / (n*n + 1))
    let p = sqrt(1.0 - z*z)
    let y = sqrt(1.0 - n*n*z*z)
    let x = sqrt(1.0 - y*y)

    print("x=\(x), y=\(y), z=\(z), p=\(p), n=\(n)")

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
        print("okokokokok1")
        let boundary = coefficients.count / 2
        print("okokokokok2")
        let zeroStateMatrix = prepareQubitStates(coefficients: Array(coefficients[0..<boundary]))
        print("okokokokok3")
        let oneStateMatrix = prepareQubitStates(coefficients: Array(coefficients[boundary..<coefficients.count]))
        print("okokokokok4")
        let firstQubit = prepareOneQubitState(
            zeroCoefficient: Array(coefficients[0..<boundary]).reduce(0, {x, y in x + pow(y, 1)}),
            oneCoefficient: Array(coefficients[boundary..<coefficients.count]).reduce(0, {x, y in x + pow(y, 1)})
        )
        print("okokokokok5")
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
        print("a1")
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
        print("a2")
        
        let result = Matrix.stackAsGates(
            [
                Matrix.kronekerProduct(matrices: [firstQubit, identity]),
                oneStateConditionedMatrix,
                zeroStateConditionedMatrix
            ]
        )
        print("a3")
        print("RESULT")
        print(result)
        print("mmmmm")
        return result
    }


    // if coefficients.count == 4 {
    //     let firstQubit = prepareOneQubitState(zeroCoefficient: pow(coefficients[0],2) + pow(coefficients[1],2), oneCoefficient: pow(coefficients[2],2) + pow(coefficients[3],2))
    //     let secondQubitL = prepareOneQubitState(zeroCoefficient: coefficients[0], oneCoefficient: coefficients[1])
    //     let secondQubitR = prepareOneQubitState(zeroCoefficient: coefficients[2], oneCoefficient: coefficients[3])
        
    //     let firstQubitInverse = Matrix.kronekerProduct(lhs: FLIP, rhs: IDENTITY)

    //     // Second qubit when first = 1
    //     let p1Part = Matrix.kronekerProduct(
    //         matrices: [
    //             P1,
    //             secondQubitR
    //         ]
    //     )
    //     let p0Part = Matrix.kronekerProduct(
    //         matrices: [
    //             P0,
    //             IDENTITY
    //         ]
    //     )
    //     let resultR = try! (
    //         p1Part + p0Part
    //     ).get()

    //     // Second qubit when first = 0
    //     let p1PartL = Matrix.kronekerProduct(
    //         matrices: [
    //             P1,
    //             secondQubitL
    //         ]
    //     )
    //     print("--")
    //     print(p1PartL)
    //     let p0PartL = Matrix.kronekerProduct(
    //         matrices: [
    //             P0,
    //             IDENTITY
    //         ]
    //     )
    //     let resultL = Matrix.multiply(matrices: [firstQubitInverse, try! (
    //         p1PartL + p0PartL
    //     ).get(), firstQubitInverse])

    //     let firstQubitResult = Matrix.kronekerProduct(matrices: [firstQubit, IDENTITY])

    //     // print(firstQubitResult.columnCount)
    //     // print(result.columnCount)
    //     // print("//")
    //     // print(firstQubitInverse)
    //     // print(resultR)
    //     // print(resultL)
    //     // print(firstQubitResult, resultR)
    //     // print(Matrix.stackAsGates([Matrix.kronekerProduct(lhs: FLIP, rhs: IDENTITY), CNOT]))
    //     return Matrix.stackAsGates([firstQubitResult, resultR, resultL])
    // }
    // return IDENTITY
}

// p1 = -0.813733471206735
// x1 = -0.8137334712067363
// y1 = 0.5812381937190964
// z1 = -0.5812381937190973

// func prepareOneQubitState(zeroCoefficient: Double, oneCoefficient: Double) -> Matrix {
//     return prepareOneQubitState(zeroCoefficient: Complex<Double>(zeroCoefficient, 0), oneCoefficient: Complex<Double>(oneCoefficient, 0))
// }
