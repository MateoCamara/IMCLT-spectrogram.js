// windows section

function windowCosine(bufferSize) {
    let windowValues = []
    for (let i = 0; i < bufferSize; i++) {
        windowValues[i] = Math.sin(Math.PI / 1024 * (i + 0.5));
    }
    return windowValues;
}

// Inverse spectrogram section

function ispectrogramMCLT(data, framelength = 1024, centered = true) {
    let overlap = 2;
    let hopsize = framelength / overlap;
    let spectrogram = data;

    let values = [];
    let chunksSize = spectrogram.length;
    for (let chunk = 0; chunk < chunksSize; chunk++) {
        values.push(math.dotMultiply(imclt(spectrogram[chunk]), windowCosine(framelength)));
    }

    let result = [];
    for (let i = 0; i < values.length; i++) {
        if (i === 0) {
            result = result.concat(values[i].slice(0, values[i].length / 2));
            result = result.concat(math.add(values[i].slice(values[i].length / 2, values[i].length), values[i + 1].slice(0, values[i + 1].length / 2)));
        } else if (i === values.length - 1) {
            result = result.concat(values[i].slice(values[i].length / 2, values[i].length));
        } else {
            result = result.concat(math.add(values[i].slice(values[i].length / 2, values[i].length), values[i + 1].slice(0, values[i + 1].length / 2)));
        }
    }

    return result;
}

// inverse transformation section

function imclt(X, odd = true) {
    if (!odd && X.length % 2 === 0) {
        throw "Even inverse CMDCT requires an odd number of coefficients"
    }

    // if (odd) {
    let N = X.length;
    let N2 = N * 2;
    let n0 = (N + 1) / 2;
    let post_twiddle = [];
    let auxMul0 = math.PI / N2;
    for (let i = 0; i < N2; i++) {
        post_twiddle.push(math.exp(math.complex(0, auxMul0 * (i + n0))));
    }

    let reversedX = math.clone(X);
    reversedX.reverse(); // reverse
    let Y1 = X;
    let Y2 = math.dotMultiply(-1, math.conj(reversedX));
    let Y = Y1.concat(Y2);
    // } else {
    //     // not odd not implemented
    // }

    let pre_twiddle = [];
    let auxMul1 = math.PI * n0 / N;
    for (let i = 0; i < N2; i++) {
        pre_twiddle.push(math.exp(math.complex(0, auxMul1 * i)));
    }

    let y = [];
    let auxMul2 = math.dotMultiply(Y, pre_twiddle);
    for (let k = 0; k < N2; k++) {
        let auxMul3 = math.PI * 2 * k / N2;
        y[k] = []
        for (let ii = 0; ii < N2; ii++) {
            y[k].push(math.dotDivide(math.dotMultiply(auxMul2[ii], math.exp(math.complex(0, auxMul3 * ii))), N2));
        }
        y[k] = math.sum(y[k]);
    }

    return math.re(math.dotMultiply(math.dotMultiply(y, post_twiddle), math.sqrt(N)));

}

// example

function example(sampleDataSize=2048) {
    let sampleData = []
    for (let i = 0; i < sampleDataSize; i++) {
        sampleData.push(i)
    }

    return ispectrogramMCLT(sampleData)
}