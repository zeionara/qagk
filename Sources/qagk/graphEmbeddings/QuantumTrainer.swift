import TensorFlow

func train(model: QuantumGraphEmbedder, lr: Float = 0.03, nEpochs: Int = 5, positiveSamples: Tensor<Int32>, negativeSamples: Tensor<Int32>) {
    for i in 1...nEpochs {
        print("Running \(i) epoch...")
        let positiveInferredLabels = model(triples: positiveSamples)
        print("Inferred label \(positiveInferredLabels)")
        let negativeInferredLabels = model(triples: negativeSamples)
        // print(derivatives)
        // print("Inferred label \(inferredLabels)")

        for layer in 0...1 {
            for qubit in 0..<model.relationshipEmbeddings[0][layer].count {
                let losses = model.computeLosses(
                    triples: positiveSamples, loss: positiveInferredLabels, layer: layer, qubit: qubit
                )
                model.updateParams(
                    triples: positiveSamples, loss: losses, lr: Double(lr), layer: layer, qubit: qubit
                )
                // model.relationshipEmbeddings[0][layer][qubit].alpha -= (trueLabel - inferredLabel) * lr * derivatives.alpha       
                // model.relationshipEmbeddings[0][layer][qubit].beta -= (trueLabel - inferredLabel) * lr * derivatives.beta
                // model.relationshipEmbeddings[0][layer][qubit].gamma -= (trueLabel - inferredLabel) * lr * derivatives.gamma
            }
        }

        let _positiveInferredLabels = model(triples: positiveSamples)
        print("Inferred label after training \(_positiveInferredLabels)")
    }
}