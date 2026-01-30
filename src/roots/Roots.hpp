#pragma once

using Function = double (*)(double);

double bisection(Function f, double a, double b, double tol = 1e-6);
double RegulaFalsi(Function f, double a, double b, double tol = 1e-6);
double newtonRaphson(Function f, Function df, double x0, double tol = 1e-6);
double secant(Function f, double x0, double x1, double tol = 1e-6);

// Finding extrema
double localmax(Function f, Function df, Function ddf, double a, double b, double tol = 1e-6);