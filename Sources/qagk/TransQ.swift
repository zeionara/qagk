import TensorFlow
import SwiftQuantumComputing

public func computeL2Norm(data: Tensor<Float>) -> Tensor<Float> {
    sqrt((data * data).sum(alongAxes: [1]))
}

public func normalizeWithL2(tensor: Tensor<Float>) -> Tensor<Float> {
    tensor / computeL2Norm(data: tensor)
}

public func initEmbeddingRotations(dimensionality: Int, nItems: Int, device device_: Device) -> Embedding<Float> {
    return Embedding(
            embeddings: normalizeWithL2(
                    tensor: Tensor<Float>(
                            randomUniform: [nItems, dimensionality],
                            lowerBound: Tensor(Float(0.0), on: device_),
                            upperBound: Tensor(Float.pi * 2, on: device_),
                            on: device_
                    )
            )
    )
}

public func rotate(x: Tensor<Float>) -> Tensor<Float> {
    let gates: [Gate] = x.unstacked().enumerated().map{(i, angle) in
        .rotation(axis: .x, radians: Double(angle.scalar!), target: i)
    }
    let circuit = MainCircuitFactory().makeCircuit(gates: gates)
    let statevector = try! circuit.statevector().get()
    print("Statevector: \(statevector)\n")
    print("Probabilities: \(statevector.probabilities())\n")
    print("Summarized probabilities: \(statevector.summarizedProbabilities())\n")
    return Tensor(0.0)
}


public struct TransQ: Module {
    public var entityEmbeddingXRotations: Embedding<Float>
    public var relationshipEmbeddingXRotations: Embedding<Float>
    // @noDerivative public let device: Device
    // @noDerivative public let dataset: KnowledgeGraphDataset<SourceElement, NormalizedElement>?

    @noDerivative public let embeddingDimensionality: Int


    public init(embeddingDimensionality: Int = 2, device device_: Device = Device.default) {
        self.embeddingDimensionality = embeddingDimensionality
        self.entityEmbeddingXRotations = initEmbeddingRotations(dimensionality: embeddingDimensionality, nItems: 3, device: device_)
        self.relationshipEmbeddingXRotations = initEmbeddingRotations(dimensionality: embeddingDimensionality, nItems: 3, device: device_)
    }

    @differentiable
    public func callAsFunction(_ triples: Tensor<Int32>) -> Tensor<Float> {
        let headEmbeddingXRotations = entityEmbeddingXRotations(triples.transposed()[0])
        for embedding in headEmbeddingXRotations.unstacked() {
            rotate(x: embedding)
        }
        // let headEmbeddings = entityEmbeddings(triples.transposed()[0])
        // print(headEmbeddings.shape)
        let tailEmbeddingXRotations = entityEmbeddingXRotations(triples.transposed()[1])
        // let tailEmbeddings = entityEmbeddings(triples.transposed()[1])
        // print(tailEmbeddings.shape)
        let relationshipEmbeddingXRotations = entityEmbeddingXRotations(triples.transposed()[2])
        // let relationshipEmbeddings_ = relationshipEmbeddings(triples.transposed()[2])
        // print(relationshipEmbeddings_.shape)
        // let score = computeScore(head: headEmbeddings, tail: tailEmbeddings, relationship: relationshipEmbeddings_)
        // print(score)
        return headEmbeddingXRotations * tailEmbeddingXRotations * relationshipEmbeddingXRotations
    }
}
