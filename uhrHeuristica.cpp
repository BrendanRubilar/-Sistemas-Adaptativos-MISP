#include <cstdint>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <sstream>

#include "utils.cpp"
#include "Metaheuristica.cpp"

int main(int argc, char *argv[])
{
    std::uint32_t base_seed = 123456789u; // base fija para reproducibilidad


    std::int64_t runs, lower, upper, step;
    validate_input(argc, argv, runs, lower, upper, step);

    int total_instances = (upper - lower) / step + 1; 
    int total_densities = 9; 
    int total_sizes = 3;     
    std::int64_t total_runs = runs * total_instances * total_densities * total_sizes;

    std::int64_t executed_runs = 0;

    // CSV global
    std::ofstream time_data(argv[1]);
    time_data << "size,density,t_mean,t_stdev,iset_mean,iset_stdev" << std::endl;

    std::cerr << "\033[0;36mRunning tests...\033[0m\n\n";

    for (int size : {1000, 2000, 3000}) { 
        std::string base_path = "/home/willi/-Sistemas-Adaptativos-MISP/dataset_grafos_no_dirigidos/new_" 
                                + std::to_string(size) + "_dataset/"; 

        // CSV por tamaño
        std::string csv_name = "resultados_" + std::to_string(size) + ".csv";
        std::ofstream size_csv(csv_name);
        size_csv << "size,density,t_mean,t_stdev,iset_mean,iset_stdev" << std::endl;

        for (int d = 1; d <= 9; ++d) {
            double density = d / 10.0;

            std::vector<double> density_times;
            std::vector<int>    density_sizes;
            density_times.reserve(total_instances * runs);
            density_sizes.reserve(total_instances * runs);

            for (int instance = lower; instance <= upper; instance += step) {
                std::ostringstream oss;
                oss << base_path << "erdos_n" << size << "_p0c0." << d
                    << "_" << instance << ".graph";
                std::string filename = oss.str();
                load_graph(filename);
                for (int r = 0; r < runs; ++r) {
                    display_progress(++executed_runs, total_runs);

                    auto begin_time = std::chrono::high_resolution_clock::now();

                    //aqui fcking aqui los test!
                    auto sol = run_grasp(filename, 10, 0.1, base_seed);

                    auto end_time = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double, std::nano> elapsed = end_time - begin_time;

                    density_times.push_back(elapsed.count());
                    density_sizes.push_back(static_cast<int>(sol.size()));
                }
            }

            // Calcular media
            auto mean_calc = [](const auto &vec)->double {
                if (vec.empty()) return 0.0;
                double s = 0.0;
                for (auto &x : vec) s += x;
                return s / vec.size();
            };

            // Calcular desviación estándar muestral
            auto stdev_calc = [](const auto &vec, double m)->double {
                size_t n = vec.size();
                if (n <= 1) return 0.0;
                double acc = 0.0;
                for (auto &x : vec) {
                    double d = x - m;
                    acc += d * d;
                }
                return std::sqrt(acc / (n - 1));
            };

            double t_mean = mean_calc(density_times);
            double iset_mean = mean_calc(density_sizes);
            double t_stdev = stdev_calc(density_times, t_mean);
            double iset_stdev = stdev_calc(density_sizes, iset_mean);

            // Escribir en CSV global
            time_data << size << ","
                      << density << ","
                      << std::setprecision(10) << t_mean << ","
                      << t_stdev << ","
                      << iset_mean << ","
                      << iset_stdev << std::endl;

            // Escribir en CSV por tamaño
            size_csv << size << ","
                     << density << ","
                     << std::setprecision(10) << t_mean << ","
                     << t_stdev << ","
                     << iset_mean << ","
                     << iset_stdev << std::endl;
        }

        size_csv.close();
    }

    std::cerr << "\n\n\033[1;32mDone!\033[0m\n";
    return 0;
}