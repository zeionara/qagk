import SwiftQuantumComputing

public extension CircuitStatevector {
    var vector: [Double] {
        var probs: [String: Double] = [String: Double]()
        for (key, value) in self.summarizedProbabilities() {
            probs[String(key.reversed())] = value
        }
        return probs.keys.sorted().map{probs[$0]!}
    }

    var firstQubitPositiveneStats: Double {
        try! self.summarizedProbabilities(byQubits: [0]).get()["1"] ?? 0
    }
}
