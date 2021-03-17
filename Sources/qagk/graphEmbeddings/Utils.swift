import Foundation

public func computeL2Norm(data: [Double]) -> Double {
    return sqrt(data.reduce(0, {x, y in x + pow(y, 2)}))
}

public func normalizeWithL2(vector: [Double]) -> [Double] {
    let norm = computeL2Norm(data: vector)
    return vector.map{$0 / norm}
}

public extension Array where Element == Double {
    var normalized: [Double] {
        normalizeWithL2(vector: self)
    }

    var asProbabilities: [Double] {
        var result = [Double]()
        for item in normalized {
            result.append(Double.pow(item, 2))
        }
        return result
    }
}
