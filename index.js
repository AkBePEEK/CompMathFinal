function AbsoluteError(x, y) {
    return Math.abs(x - y)
}

function BisectionMethod(f, a, b, tol){
    if (f(a) * f(b) >= 0){
        console.log("Invalid initial values. f(a) and f(b) must be of different signs.")
        return null
    }

    let midpoint = (a + b) / 2

    while (Math.abs(f(midpoint)) > tol){
        if (f(a) * f(midpoint) < 0)
            b = midpoint
        else {
            a = midpoint
            midpoint = (a + b) / 2
        }
    }
    return midpoint
}

function SecantMethod(f, x0, x1, tol){
    while (Math.abs(f(x1)) > tol){
        let x_temp = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        x0 = x1
        x1 = x_temp
    }
    return x1
}