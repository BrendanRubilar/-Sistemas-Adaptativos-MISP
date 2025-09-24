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

// ESTA SHIT ESTA COMENTADA PORQUE GENERA EL CSV DE ESTO COMPLETO: "/home/franc/-Sistemas-Adaptativos-MISP/dataset_grafos_no_dirigidos/new_1000_dataset/"; 
// HAY QUE VERIFICAR DE ALGUNA FORMA? QUE LOS RESULTADOS SON CORRECTOS. PERO EL MAIN FUNCIONAL ES EL QUE DEMUESTRA QUE LEE CORRECTAMENTE LOS ARCHVIOS DE LOS GRAFOS C:

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

            std::ostringstream oss;
            std::string base_path = "/home/willi/-Sistemas-Adaptativos-MISP/dataset_grafos_no_dirigidos/new_3000_dataset/"; 
            oss << base_path << "erdos_n3000_p0c0.";

            // Manejar los diferentes formatos de densidad de forma individual para evitar errores de precisión por instancias como 0.05, 0.15, 0.25, etc.
            if (std::abs(density - 0.05) < 1e-9) {
                oss << "05";
            } else if (std::abs(density - 0.1) < 1e-9) {
                oss << "1";
            } else if (std::abs(density - 0.15) < 1e-9) {
                oss << "15";
            } else if (std::abs(density - 0.2) < 1e-9) {
                oss << "2";
            } else if (std::abs(density - 0.25) < 1e-9) {
                oss << "25";
            } else if (std::abs(density - 0.3) < 1e-9) {
                oss << "3";
            } else if (std::abs(density - 0.35) < 1e-9) {
                oss << "35";
            } else if (std::abs(density - 0.4) < 1e-9) {
                oss << "4";
            } else if (std::abs(density - 0.45) < 1e-9) {
                oss << "45";
            } else if (std::abs(density - 0.5) < 1e-9) {
                oss << "5";
            } else if (std::abs(density - 0.55) < 1e-9) {
                oss << "55";
            } else if (std::abs(density - 0.6) < 1e-9) {
                oss << "6";
            } else if (std::abs(density - 0.65) < 1e-9) {
                oss << "65";
            } else if (std::abs(density - 0.7) < 1e-9) {
                oss << "7";
            } else if (std::abs(density - 0.75) < 1e-9) {
                oss << "75";
            } else if (std::abs(density - 0.8) < 1e-9) {
                oss << "8";
            } else if (std::abs(density - 0.85) < 1e-9) {
                oss << "85";
            } else if (std::abs(density - 0.9) < 1e-9) {
                oss << "9";
            } else if (std::abs(density - 0.95) < 1e-9) {
                oss << "95";
            }

            oss << "_" << instance << ".graph";
            std::string filename = oss.str();
            // Múltiples ejecuciones para estadística
            for (i = 0; i < runs; i++) {
                display_progress(++executed_runs, total_runs);

                begin_time = std::chrono::high_resolution_clock::now();
                // ==== FUNCIONES A TESTEAR ====
                
                load_graph(filename);
                //auto iset = misp_heuristica_greedy();
                auto iset = misp_greedy_aleatorizado(0.1);
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