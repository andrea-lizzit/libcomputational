#include "utilities.h"

float* linspace(float a, float b, int n)
{
    float* arr = new float[n + 1];
    for (int i = 0; i <= n; ++i) {
        arr[i] = a + (b - a) / n * i;
    }
    return arr;
}
double* linspace(double a, double b, int n)
{
    double* arr = new double[n + 1];
    for (int i = 0; i <= n; ++i) {
        arr[i] = a + (b - a) / n * i;
    }
    return arr;
}

float* range(float offset, float h, int n)
{
    float* arr = new float[n];
    for (int i = 0; i < n; ++i) {
        arr[i] = offset + h * i;
    }
    return arr;
}

double* range(double offset, double h, int n)
{
    double* arr = new double[n];
    for (int i = 0; i < n; ++i) {
        arr[i] = offset + h * i;
    }
    return arr;
}