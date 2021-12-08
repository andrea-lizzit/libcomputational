#ifndef INTEGRALS
#define INTEGRALS
#include <functional>
#include "utilities.h"

double integrate_naif(double a, double b, int n, std::function<double(double)> fun);
double integrate_rect(double a, double b, int n, std::function<double(double)> fun);
double integrate_trap(double a, double b, int n, std::function<double(double)> fun);
double integrate_simpson(double a, double b, int n, std::function<double(double)> fun);


#endif