#include "Roots.hpp"
#include <cmath>
#include <stdexcept>

// Implements the bisection method to find a root of the function f in the interval [a, b] with a specified tolerance tol
double bisection(Function f, double a, double b, double tol)
{
    // evaluate function at endpoints
    double fa = f(a);
    double fb = f(b);

    // Check that f(a) and f(b) have opposite signs as bisection method requires opposite signs
    if (fa * fb > 0)
        throw std::runtime_error("f(a) and f(b) must have opposite signs");

    // Main loop, iterating 1 mil times at max if error
        for (int i= 0; i<1e6; ++i){
            // Locates and evaluates at midpoint
            double c = 0.5 *(a+b);
            double fc = f(c);

            // checks if the function value at midpoint is near a root
            if (std::abs(fc) < tol || (b-a)/2 <tol)
                return c; // root found

            // check if root is in left or right subinterval
            if (fa * fc < 0){
                b = c;
                fb = fc;
            }
            else {
                a = c;
                fa = fc;
            }
        }
        throw std::runtime_error("Bisection method did not converge within the maximum number of iterations");
}

// Implements the Regula Falsi (False Position) method to find a root of the function f in the interval [a, b] with a specified tolerance tol
double regula_falsi(Function f, double a, double b, double tol)
{
    // evaluate function at endpoints
    double fa = f(a);
    double fb = f(b);

    // Check that f(a) and f(b) have opposite signs as Regula Falsi method requires opposite signs
    if (fa * fb > 0)
        throw std::runtime_error("f(a) and f(b) must have opposite signs");

    // Main loop, iterating 1 mil times at max if error
    for (int i = 0; i < 1e6; ++i) {
        // Calculate the point c using the Regula Falsi formula
        double c = (a * fb - b * fa) / (fb - fa);
        double fc = f(c);

        // Find if the function value at c is near a root
        if (std::abs(fc) < tol)
            return c; // root found

        // Update the interval [a, b] based on the sign of f(c)
        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    throw std::runtime_error("Regula Falsi method did not converge within the maximum number of iterations");
}

// Implements the Newton-Raphson method to find a root of the function f given its derivative df, starting from an initial guess x0 with a specified tolerance tol
double newton_raphson(Function f, Function df, double x0, double tol)
{

    // Initial guess
    double x = x0;

    // Main loop, iterating 1 mil times at max if error
    for (int i = 0; i < 1e6; ++i) {
        double fx = f(x);
        double dfx = df(x);

        // Check if derivative is zero to avoid division by zero
        if (std::abs(dfx) <= 1e-12)
            throw std::runtime_error("Derivative is zero. No solution found.");

        // Use Newton-Raphson formula to find new xn
        double xn = x - fx / dfx;

        // Check for convergence, aka xn - x < 1e-6
        if (std::abs(xn - x) < tol)
            return xn; // root found

        x = xn;
    }
    throw std::runtime_error("Newton-Raphson method did not converge within the maximum number of iterations");
}

//Implements the Secant method to find a root of the function f, starting from two initial guesses x0 and x1 with a specified tolerance tol
double secant(Function f, double x0, double x1, double tol)
{
    // Define function values at initial guesses
    double f0 = f(x0);
    double f1 = f(x1);

    for (int i =0; i < 1e6; ++i){
        // Check if division is invalid
        if (std::abs(f1 - f0) <= 1e-12)
            throw std::runtime_error("Division by zero in Secant method. No solution found.");

        // Secant method formula to find new xn1
        double xn1 = x0 - f0 * (x1 - x0) / (f1 - f0);

        //Check for convergence
        if (std::abs(xn1 - x1) < tol)
            return xn1; // root found


        // Update variables for next iteration
        x0 = x1;
        f0 = f1;
        x1 = xn1;
        f1 = f(x1);
    }
    throw std::runtime_error("Secant method did not converge within the maximum number of iterations");
}


// Find local max using where f'(x) = 0, confirm with  f''(x) <0
double localmax(Function f, Function df, Function ddf, double a, double b, double tol)  {
    double x = bisection(df, a, b, tol); // find critical point where f'(x) = 0

    if (ddf(x) < 0)
        return x; // local maximum found
    else
        throw std::runtime_error("No local maximum found in the given interval");
}
