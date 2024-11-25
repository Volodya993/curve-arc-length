#include <iostream>
#include <cmath>
#include <functional>

double integrate(std::function<double(double)> f, double a, double b, double epsilon) {
    double n = 1; 
    double h = (b - a) / n;
    double integral = 0.0;

    for (int i = 0; i < n; ++i) {
        integral += (f(a + i * h) + f(a + (i + 1) * h)) * h / 2.0;
    }

    double prev_integral;
    do {
        prev_integral = integral;
        n *= 2; 
        h = (b - a) / n;
        integral = 0.0;

        for (int i = 0; i < n; ++i) {
            integral += (f(a + i * h) + f(a + (i + 1) * h)) * h / 2.0;
        }
    } while (std::abs(integral - prev_integral) > epsilon);

    return integral;
}

double derivative(double x) {
    return -1 / x; 
}

double arcLengthFunction(double x) {
    double y_prime = derivative(x);
    return std::sqrt(1 + y_prime * y_prime);
}

int main() {
    double a = std::sqrt(3);
    double b = std::sqrt(8);
    double epsilon = 1e-6;

    double length = integrate(arcLengthFunction, a, b, epsilon);

    std::cout << "length of the curve arc is " << length << std::endl;

    return 0;
}
