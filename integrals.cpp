#include <functional>
#include <assert.h>
#include "utilities.h"

double integrate_naif(double a, double b, int n, std::function<double(double)> fun)
{
    double* domain = linspace(a, b, n);
    double integral = 0;
    double h = (b - a) / n;
    for (int i = 0; i < n; ++i) {
        integral += fun(domain[i]) * h;
    }
    free(domain);
    return integral;
}

double integrate_rect(double a, double b, int n, std::function<double(double)> fun)
{
    double* domain = linspace(a, b, n);
    double integral = 0;
    double h = (b - a) / n;
    for (int i = 0; i < n; ++i) {
        double midpoint = (domain[i] + domain[i+1]) / 2;
        integral += fun(midpoint) * h;
    }
    free(domain);
    return integral;
}
double integrate_trap(double a, double b, int n, std::function<double(double)> fun)
{
    // make a grid with n+2 points; we need n+1 points to apply the trapezoid method;
    // being a closed method, the last point corresponds to infinity and can't be evaluated (produces Nan)
    // therefore we make and n+2 grid and ignore the last value
    // the approximation is good for high n values
    // but goes asymptotically as 1/N, like naif rectangles
    double* domain = linspace(a, b, n+1);
    double integral = 0;
    double h = (b - a) / n;
    integral += (fun(domain[0]) + fun(domain[n])) * h / 2;
    for (int i = 1; i < n; ++i) {
        integral += fun(domain[i]) * h;
    }
    free(domain);
    return integral;
}

double integrate_simpson(double a, double b, int n, std::function<double(double)> fun)
{
    assert(n % 2 == 0);
    // make a grid with n+2 points as for the trapezoid method
    double* domain = linspace(a, b, n+1);
    double integral = 0;
    double h = (b - a) / n;
    integral += (fun(domain[0]) + fun(domain[n])) * h / 3;
    for (int i = 1; i < n; ++i) {
        integral += fun(domain[i]) * h / 3 * (i % 2 ? 4 : 2);
    }
    free(domain);
    return integral;
}