#include <iostream>
#include <vector>
#include <cmath>
#include <string>
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

std::vector<double> finiteDifferenceNeumann(double a, double b, double ua_d, double ub, int n, std::function<double(double)> f, int method) {
    double h = (b - a) / (n + 1);
    int m = n + 1;

    std::vector<double> diag(m, -2.0 / pow(h, 2));
    std::vector<double> lower(m - 1, 1.0 / pow(h, 2));
    std::vector<double> upper(m - 1, 1.0 / pow(h, 2));
    diag[0] = -1.0 / h;
    upper[0] = 1.0 / h;

    std::vector<double> F(m);
    for (int i = 0; i < m; i++) {
        double x = a + (i + 1) * h;
        F[i] = f(x);
    }
    F[0] = ua_d;

    if (method == 2) {
        F[0] += (h / 2) * f(a);
    }

    F[m - 1] -= ub / pow(h, 2);

    std::vector<double> y(m);
    thomasAlgorithm(lower, diag, upper, F, y);

    std::vector<double> r(n + 2);
    r[n + 1] = ub;
    for (int i = 0; i < m; i++) {
        r[i] = y[i];
    }

    return r;
}

double pNorm(const std::vector<double>& vec, int p) {
    double sum = 0.0;
    for (double v : vec) {
        sum += std::pow(v, p);
    }
    return std::pow(sum, 1.0 / p);
}

double calculateEOC(const std::vector<double>& errorCoarse, const std::vector<double>& errorFine, double p) {
    double numerator = pNorm(errorCoarse, p);
    double denominator = pNorm(errorFine, p);
    return std::log2(numerator / denominator);
}

int main() {
    double a = 0;
    double b = 1;
    double ua_d = 1; // Neumannova OP
    double ub = 4;  // Dirichletova OP
    std::function<double(double)> f = [](double t) { return 12 * t * t + 2; };

    std::function<double(double)> analyticalSolution = [](double t) { return t * t * t * t + t * t + t + 1; };

    int method = 2;
    std::vector<int> gridSizes = { 10, 20, 40, 80, 160 };

    std::vector<double> previousSolution;
    std::vector<double> previousGridPoints;
    double previousH = 0;

    for (int i = 0; i < gridSizes.size(); i++) {
        int n = gridSizes[i];
        double h = (b - a) / (n + 1);
        std::vector<double> currentSolution = finiteDifferenceNeumann(a, b, ua_d, ub, n, f, method);
        std::vector<double> currentGridPoints(n + 2);
        for (int j = 0; j <= n + 1; j++) {
            currentGridPoints[j] = a + j * h;
        }

        std::vector<double> currentError(n + 2);
        for (int j = 0; j <= n + 1; j++) {
            double t = currentGridPoints[j];
            currentError[j] = std::abs(currentSolution[j] - analyticalSolution(t));
        }

        if (i > 0) {
            std::vector<double> interpolatedErrors(previousGridPoints.size());
            for (size_t k = 0; k < previousGridPoints.size(); k++) {
                auto closest = std::min_element(currentGridPoints.begin(), currentGridPoints.end(),
                    [p = previousGridPoints[k]](double a, double b) { return std::abs(a - p) < std::abs(b - p); });
                size_t index = std::distance(currentGridPoints.begin(), closest);
                interpolatedErrors[k] = currentError[index];
            }
            double eoc = calculateEOC(interpolatedErrors, currentError, 2);
            std::cout << "EOC: " << eoc << std::endl;
        }

        std::cout << "Grid Size: " << n << std::endl;
        for (int j = 0; j <= n + 1; j++) {
            double t = a + j * h;
            double exact = analyticalSolution(t);
            double num = currentSolution[j];
            double localError = std::abs(num - exact);
            std::cout << "t = " << t << ", Numerical = " << num << ", Exact = " << exact << ", Error = " << localError << std::endl;
        }

        previousSolution = currentSolution;
        previousGridPoints = currentGridPoints;
        previousH = h;
    }

    return 0;
}