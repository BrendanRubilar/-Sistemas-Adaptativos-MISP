#include <cstdint>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <iomanip> // para controlar formato decimal

#include "utils.cpp"
#include "Heuristica.cpp"

int main(int argc, char *argv[])
{
    // Validate and sanitize input
    std::int64_t runs, lower, upper, step;
    validate_input(argc, argv, runs, lower, upper, step);

    // Clock vars
    std::int64_t i, executed_runs;
    // Total de corridas (densidades * instancias * runs)
    int total_instances = (upper - lower) / step + 1; // si lower=1, upper=30, step=1 → 30
    int total_densities = 9; // de 0.1 a 0.9
    std::int64_t total_runs = runs * total_instances * total_densities;

    std::vector<double> times(runs);
    double mean_time, time_stdev, dev;
    auto begin_time = std::chrono::high_resolution_clock::now();
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::nano> elapsed_time = end_time - begin_time;

    // RNG
    std::random_device rd;
    std::mt19937_64 rng(rd());

    // Output file
    std::ofstream time_data(argv[1]);
    time_data << "density,instance,t_mean,t_stdev" << std::endl;

    // Run tests
    std::cerr << "\033[0;36mRunning tests...\033[0m" << std::endl << std::endl;
    executed_runs = 0;

    for (int d = 1; d <= 9; d++) { // densidades 0.1 ... 0.9
        double density = d / 10.0;

        for (int instance = lower; instance <= upper; instance += step) {
            mean_time = 0;
            time_stdev = 0;

            // Nombre del archivo: erdos_n1000_p0c0.X_Y.graph
            std::ostringstream oss;
            oss << "erdos_n1000_p0c0_" << std::fixed << std::setprecision(1)
                << density << "_" << instance << ".graph";
            std::string filename = oss.str();

            // Múltiples ejecuciones para estadística
            for (i = 0; i < runs; i++) {
                display_progress(++executed_runs, total_runs);

                begin_time = std::chrono::high_resolution_clock::now();
                // ==== FUNCIONES A TESTEAR ====
                load_graph_limit_nodes(filename, 1000);
                auto orden = nodos_ordenados_por_grado_simple();
                auto iset = misp_heuristica_por_orden(orden);
                // ==============================
                end_time = std::chrono::high_resolution_clock::now();

                elapsed_time = end_time - begin_time;
                times[i] = elapsed_time.count();

                mean_time += times[i];
            }

            // Media
            mean_time /= runs;

            // Desviación estándar
            for (i = 0; i < runs; i++) {
                dev = times[i] - mean_time;
                time_stdev += dev * dev;
            }
            time_stdev = std::sqrt(time_stdev / (runs - 1));

            // Escribir fila con resultados
            time_data << density << "," << instance << ","
                      << mean_time << "," << time_stdev << std::endl;
        }
    }

    std::cerr << std::endl << std::endl;
    std::cerr << "\033[1;32mDone!\033[0m" << std::endl;

    time_data.close();
    return 0;
}
