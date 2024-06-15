#include <iostream>
#include <cmath>
#include <limits>

using namespace std;

// Function to define the tolerance for convergence
double tolerance() {
    return numeric_limits<double>::epsilon() * 100; // Adjust based on desired accuracy
}

// Bisection Method implementation
class BisectionMethod {
public:
    BisectionMethod(double a, double b) : lowerBound(a), upperBound(b) {
        validateInput();
    }

    double findRoot() {
        double root, fa, fb;

        while (abs(upperBound - lowerBound) >= tolerance()) {
            root = (lowerBound + upperBound) / 2.0;
            fa = evaluateFunction(lowerBound);
            fb = evaluateFunction(upperBound);

            if (fa * fb < 0) { // Root lies in the interval
                upperBound = root;
            } else {
                lowerBound = root;
            }
        }

        return root;
    }

private:
    double lowerBound, upperBound;

    // Replace with your actual function definition here
    double evaluateFunction(double x) {
        // Example function: You can replace this with your actual function
        return x * x * x - x - 1; // x^3 - x - 1
    }

    void validateInput() {
        if (evaluateFunction(lowerBound) * evaluateFunction(upperBound) >= 0) {
            throw invalid_argument("Function values at lower and upper bounds must have opposite signs.");
        }
    }
};

// False Position Method implementation
class FalsePositionMethod {
public:
    FalsePositionMethod(double a, double b) : lowerBound(a), upperBound(b) {
        validateInput();
    }

    double findRoot() {
        double root, fa, fb;

        while (abs(upperBound - lowerBound) >= tolerance()) {
            fa = evaluateFunction(lowerBound);
            fb = evaluateFunction(upperBound);

            // Calculate the false position (root estimate)
            root = lowerBound - (fa * (upperBound - lowerBound)) / (fb - fa);

            if (abs(evaluateFunction(root)) <= tolerance()) { // Root found with sufficient accuracy
                return root;
            } else if (fa * evaluateFunction(root) < 0) {
                upperBound = root;
            } else {
                lowerBound = root;
            }
        }

        return root; // Approximate root within tolerance
    }

private:
    double lowerBound, upperBound;

    // Replace with your actual function definition here
    double evaluateFunction(double x) {
        // Example function: You can replace this with your actual function
        return x * x * x - 4 * x + 1; // x^3 - 4x + 1
    }

    void validateInput() {
        if (evaluateFunction(lowerBound) * evaluateFunction(upperBound) >= 0) {
            throw invalid_argument("Function values at lower and upper bounds must have opposite signs.");
        }
    }
};

// Newton-Raphson Method implementation
class NewtonRaphsonMethod {
public:
    NewtonRaphsonMethod(double initialGuess) : x0(initialGuess) {}

    double findRoot() {
        double currentGuess, previousGuess, functionValue, derivativeValue;

        currentGuess = x0;
        previousGuess = 0; // Initialize for first iteration

        while (abs(currentGuess - previousGuess) >= tolerance() * abs(currentGuess)) {
            previousGuess = currentGuess;
            functionValue = evaluateFunction(currentGuess);
            derivativeValue = differentiateFunction(currentGuess);

            // Check for division by zero (derivative is very close to zero)
            if (abs(derivativeValue) < tolerance()) {
                throw domain_error("Derivative is close to zero, possibly causing division by zero.");
            }

            // Update current guess using the Newton-Raphson formula
            currentGuess = currentGuess - functionValue / derivativeValue;
        }

        return currentGuess;
    }

private:
    double x0;

    // Replace with your actual function definition here
    double evaluateFunction(double x) {
        // Example function: You can replace this with your actual function
        return x - cos(x); // x - cos(x)
    }

    //
// Replace with your actual derivative definition here
double differentiateFunction(double x) {
    return 1 + sin(x); // Example derivative: 1 + sin(x)
}
};

// Secant Method implementation
class SecantMethod {
public:
    SecantMethod(double initialGuess1, double initialGuess2) : x0(initialGuess1), x1(initialGuess2) {
        validateInput();
    }

    double findRoot() {
        double currentGuess, previousGuess;

        currentGuess = x1;
        previousGuess = x0;

        while (abs(currentGuess - previousGuess) >= tolerance()) {
            // Calculate the secant formula for the next approximation
            currentGuess = currentGuess - evaluateFunction(currentGuess) * (currentGuess - previousGuess) / (evaluateFunction(currentGuess) - evaluateFunction(previousGuess));
            previousGuess = currentGuess;
        }

        return currentGuess;
    }

private:
    double x0, x1;

    // Replace with your actual function definition here
    double evaluateFunction(double x) {
        // Example function: You can replace this with your actual function
        return x * x * x + 4 * x * x - 10; // x^3 + 4x^2 - 10
    }

    void validateInput() {
        if (abs(x0 - x1) <= tolerance()) {
            throw invalid_argument("Initial guesses must be distinct for the Secant Method.");
        }
    }
};

// Fixed Point Method implementation
class FixedPointMethod {
public:
    FixedPointMethod(double initialGuess) : x0(initialGuess) {}

    double findRoot() {
        double currentGuess, previousGuess, absoluteError;

        currentGuess = x0;
        previousGuess = 0; // Initialize for first iteration

        while (abs(currentGuess - previousGuess) >= tolerance()) {
            previousGuess = currentGuess;
            currentGuess = iterateFunction(previousGuess);

            absoluteError = abs(currentGuess - previousGuess);
        }

        return currentGuess;
    }

private:
    double x0;

    // Replace with your actual function definition for g(x)
    double iterateFunction(double x) {
        // Example function: You can replace this with your actual function definition
        return exp(-x); // e^-x
    }
};

int main() {
    // Example usage:
    try {
        // Bisection Method
        double a = 1, b = 2;
        BisectionMethod bm(a, b);
        double root = bm.findRoot();
        cout << "Bisection Method: The approximate root is: " << root << endl;

        // False Position Method
        a = 1, b = 2;
        FalsePositionMethod fm(a, b);
        root = fm.findRoot();
        cout << "False Position Method: The approximate root is: " << root << endl;

        // Newton-Raphson Method
        double x0 = 1;
        NewtonRaphsonMethod nr(x0);
        root = nr.findRoot();
        cout << "Newton-Raphson Method: The approximate root is: " << root << endl;

        // Secant Method
        double x1 = 2;
        SecantMethod sm(x0, x1);
        root = sm.findRoot();
        cout << "Secant Method: The approximate root is: " << root << endl;

        // Fixed Point Method
        double A0 = 1;
        FixedPointMethod fpi(A0);
        root = fpi.findRoot();
        cout << "Fixed Point Method: The approximate root is: " << root << endl;
    } catch (const invalid_argument& e) {
        cerr << "Error: " << e.what() << endl;
    } catch (const domain_error& e) {
        cerr << "Error: " << e.what() << endl;
    }

    return 0;
}
