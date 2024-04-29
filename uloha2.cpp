#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>

double pNorm(const std::vector<double>& vec, double h, int p) {
    double sum = 0.0;
    for (double v : vec) {
        sum += std::pow(h * std::fabs(v), p);
    }
    return std::pow(sum, 1.0 / p);
}

double calculateEOC(const std::vector<double>& errorCoarse, const std::vector<double>& errorFine, double h, int p) {
    double numerator = pNorm(errorCoarse, h, p);
    double denominator = pNorm(errorFine, h / 2, p);
    return std::log2(numerator / denominator);
}

std::vector<double> solve_fdm(double h, int N) {
    std::vector<double> u(N + 1, 0.0);
    u[0] = 0;
    u[1] = u[0] + h;

    for (int i = 1; i < N; ++i) {
        double t_i = i * h;
        if (i < N - 1) {
            u[i + 1] = 2 * u[i] - u[i - 1] + h * h * (12 * t_i * t_i + 2);
        }
    }
    u[N] = 4;

    return u;
}

double exact_solution(double t) {
    return 1 + t + pow(t, 2) + pow(t, 4);
}

void write_to_csv(const std::string& filename, std::vector<double>& u, double h, int N) {
    std::ofstream file(filename);
    file << "t,u_numeric,u_exact,error\n";
    for (int i = 0; i <= N; ++i) {
        double t_i = i * h;
        if (i != N) u[i]++;
        file << t_i << "," << u[i] << "," << exact_solution(t_i) << "," << fabs(u[i] - exact_solution(t_i)) << "\n";
    }
    file.close();
}

int main() {
    std::vector<int> hs = { 5, 10, 20, 40, 80 };
    std::vector<double> previous_solution, current_errors, previous_errors;
    std::vector<std::pair<double, double>> eoc_results;

    for (int N : hs) {
        double h = 1.0 / N;
        std::vector<double> current_solution = solve_fdm(h, N);
        current_errors.clear();

        for (int i = 0; i <= N; ++i) {
            current_errors.push_back(std::fabs(current_solution[i] - exact_solution(i * h)));
        }

        write_to_csv("solution_h_" + std::to_string(N) + ".csv", current_solution, h, N);

        if (!previous_errors.empty()) {
            double EOC = calculateEOC(previous_errors, current_errors, h, 2);
            eoc_results.emplace_back(h, EOC);
            std::cout << "EOC for h = " << h << ": " << EOC << std::endl;
        }

        previous_solution = current_solution;
        previous_errors = current_errors;
    }

    std::ofstream eoc_file("eoc_results.csv");
    eoc_file << "h,EOC\n";
    for (const auto& result : eoc_results) {
        eoc_file << result.first << "," << result.second << "\n";
    }
    eoc_file.close();

    return 0;
}
