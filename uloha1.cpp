#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>
#include <functional>

#define Pi 3.14159265358979323846

void thomasAlgorithm(const std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d, std::vector<double>& u) {
    int n = d.size() - 1;

    for (int i = 1; i < n; i++) {
        double m = a[i] / b[i - 1];
        b[i] = b[i] - m * c[i - 1];
        d[i] = d[i] - m * d[i - 1];
    }

    u[n] = d[n] / b[n];
    for (int i = n - 1; i >= 0; i--) {
        u[i] = (d[i] - c[i] * u[i + 1]) / b[i];
    }
}

double pNorm(const std::vector<double>& vec, int p) {
    double sum = 0.0;
    for (double v : vec) {
        sum += std::pow(v, p);
    }
    return std::pow(sum, 1.0 / p);
}

void solveFiniteDifferences(int n, std::vector<double>& u, double& h) {
    h = Pi / (n - 1);
    std::vector<double> a(n, 1.0), b(n, -2.0), c(n, 1.0), d(n, 0.0);

    for (int i = 1; i < n - 1; i++) {
        d[i] = -h * h * std::sin(i * h);
    }

    d[0] = 1.0;
    d[n - 1] = 1.0;
    a[0] = 0.0;
    c[n - 1] = 0.0;

    thomasAlgorithm(a, b, c, d, u);
}

double calculateEOC(const std::vector<double>& errorCoarse, const std::vector<double>& errorFine, double p) {
    double numerator = pNorm(errorCoarse, p);
    double denominator = pNorm(errorFine, p);
    return std::log2(numerator / denominator);
}

int main() {
    std::vector<int> steps = { 10, 20, 40, 80, 160 };

    std::function<double(double)> exactSolutionFunction = [](double t) -> double { return 1 + std::sin(t); };

    std::vector<double> previousSolution;
    double previousH = Pi;
    double eocP1 = 0.0, eocP2 = 0.0;

    for (int i = 0; i < steps.size(); i++) {
        int n = steps[i];
        double h = Pi / (n - 1);
        std::vector<double> u(n), error(n);

        solveFiniteDifferences(n, u, h);

        for (int j = 0; j < n; j++) {
            double exact = exactSolutionFunction(j * h);
            error[j] = std::abs(u[j] - exact);
        }

        if (i > 0) {
            double normErrorP1 = h * pNorm(error, 1);
            double normErrorP2 = std::sqrt(h) * pNorm(error, 2);

            eocP1 = std::log2(previousH / h) * std::log2(pNorm(previousSolution, 1) / normErrorP1);
            eocP2 = std::log2(previousH / h) * std::log2(pNorm(previousSolution, 2) / normErrorP2);

            std::cout << "EOC for p=1: " << eocP1 << std::endl;
            std::cout << "EOC for p=2: " << eocP2 << std::endl;
        }

        previousSolution = error;
        previousH = h;
    }

    return 0;
}
