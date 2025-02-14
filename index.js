function absoluteError(x, y) {
    return Math.abs(x - y);
}

function relativeError(x, y) {
    return Math.abs(x - y) / y;
}

function bisectionMethod(f, a, b, tol){
    if (f(a) * f(b) >= 0){
        console.log("Invalid initial values. f(a) and f(b) must be of different signs.");
        return null;
    }

    let midpoint = (a + b) / 2;
    let iteration = 0
    while (Math.abs(f(midpoint)) > tol){
        if (f(a) * f(midpoint) < 0)
            b = midpoint;
        else {
            a = midpoint;
            midpoint = (a + b) / 2;
        }
        iteration += 1;
    }
    return [midpoint, iteration];
}

function secantMethod(f, x0, x1, tol){
    let iteration = 0;
    while (Math.abs(f(x1)) > tol){
        let x_temp = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
        x0 = x1;
        x1 = x_temp;
        iteration += 1;
    }
    return [x1, iteration];
}

function jacobiMethod(coefficients, initialGuess, maxIterations = 100, tolerance = 1e-6) {
    let [x, y, z] = initialGuess;

    for (let iteration = 0; iteration < maxIterations; iteration++) {
        // Save previous values
        let xPrev = x, yPrev = y, zPrev = z;

        // Update values using the Jacobi Method
        x = (coefficients[0].constant - coefficients[0].y * yPrev - coefficients[0].z * zPrev) / coefficients[0].x;
        y = (coefficients[1].constant - coefficients[1].x * xPrev - coefficients[1].z * zPrev) / coefficients[1].y;
        z = (coefficients[2].constant - coefficients[2].x * xPrev - coefficients[2].y * yPrev) / coefficients[2].z;

        // Check for convergence
        if (
            Math.abs(x - xPrev) < tolerance &&
            Math.abs(y - yPrev) < tolerance &&
            Math.abs(z - zPrev) < tolerance
        ) {
            return { x, y, z, iterations: iteration + 1 };
        }
    }

    throw new Error("The method did not converge within the maximum number of iterations.");
}

function iterativeInverse(A, B, tol, maxIter) {
    const n = A.length;
    const I = Array.from({ length: n }, (_, i) =>
        Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))
    ); // Create the identity matrix

    for (let iter = 0; iter < maxIter; iter++) {
        // Calculate the error matrix E = A * B - I
        const E = matrixMultiply(A, B);
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                E[i][j] -= I[i][j];
            }
        }

        // Update B: B_new = B - B * E
        const B_new = matrixSubtract(B, matrixMultiply(B, E));

        // Check for convergence using the Frobenius norm of E
        let frobeniusNorm = 0;
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                frobeniusNorm += E[i][j] ** 2;
            }
        }
        frobeniusNorm = Math.sqrt(frobeniusNorm);

        if (frobeniusNorm < tol) {
            return B_new; // Return the inverse if convergence is achieved
        }

        B = B_new; // Update B for the next iteration
    }

    return B; // Return the result after maxIter iterations
}

// Helper function to multiply two matrices
function matrixMultiply(A, B) {
    const n = A.length;
    const result = Array.from({ length: n }, () => Array(n).fill(0));
    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
            for (let k = 0; k < n; k++) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return result;
}

// Helper function to subtract two matrices
function matrixSubtract(A, B) {
    const n = A.length;
    const result = Array.from({ length: n }, () => Array(n).fill(0));
    for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
            result[i][j] = A[i][j] - B[i][j];
        }
    }
    return result;
}

function plotLinearFit(x, y) {
    // Fit a straight line (linear regression)
    const xMean = math.mean(x);
    const yMean = math.mean(y);

    // Calculate slope (m) and intercept (c)
    const numerator = math.sum(math.dotMultiply(math.subtract(x, xMean), math.subtract(y, yMean)));
    const denominator = math.sum(math.dotMultiply(math.subtract(x, xMean), math.subtract(x, xMean)));
    const m = numerator / denominator; // Slope
    const c = yMean - m * xMean; // Intercept

    // Line equation: y = mx + c
    const fitLine = math.add(math.multiply(x, m), c);

    // Plot data points and best-fit line
    const trace1 = {
        x: x,
        y: y,
        mode: 'markers',
        type: 'scatter',
        name: 'Data Points',
        marker: { color: 'red' }
    };

    const trace2 = {
        x: x,
        y: fitLine,
        mode: 'lines',
        type: 'scatter',
        name: `Best Fit: y = ${m.toFixed(2)}x + ${c.toFixed(2)}`,
        line: { color: 'blue' }
    };

    const layout = {
        title: 'Graphical Method: Linear Fit',
        xaxis: { title: 'X-axis' },
        yaxis: { title: 'Y-axis' },
        showlegend: true
    };

    Plotly.newPlot('plot', [trace1, trace2], layout);
}

function newtonsForwardDifference(x, y, x_target) {
    let n = y.length;
    let h = x[1] - x[0];

    // Create forward difference table
    let forwardDiff = Array(n).fill(0).map(() => Array(n).fill(0));

    for (let i = 0; i < n; i++) {
        forwardDiff[i][0] = y[i];
    }

    for (let j = 1; j < n; j++) {
        for (let i = 0; i < n - j; i++) {
            forwardDiff[i][j] = forwardDiff[i + 1][j - 1] - forwardDiff[i][j - 1];
        }
    }

    // First derivative using Newton's forward formula at x = x_target
    let index = x.indexOf(x_target);
    if (index === -1) {
        console.log("x_target not found in x array");
        return;
    }

    let dy_dx = forwardDiff[index][1] / h;

    console.log(`First derivative at x = ${x_target}: ${dy_dx}`);
}