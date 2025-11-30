#include <cstdint>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <sstream>

#include "utils.cpp"
#include "Heuristica.cpp"

int main(int argc, char *argv[])
{
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

                for (int r = 0; r < runs; ++r) {
                    display_progress(++executed_runs, total_runs);

                    auto begin_time = std::chrono::high_resolution_clock::now();

                    load_graph(filename);
                    //auto iset = misp_heuristica_greedy();
                    auto iset = misp_greedy_aleatorizado(0.1f);

                    auto end_time = std::chrono::high_resolution_clock::now();
                    std::chrono::duration<double, std::nano> elapsed = end_time - begin_time;

                    density_times.push_back(elapsed.count());
                    density_sizes.push_back((int)iset.size());
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

*/



// Main de lectura de grafos, el funcionamiento esta arriba, LEE MUY LENTO LOS GRAFOS DE 1000 NODOS.
namespace fs = std::filesystem;

int main(int argc, char *argv[]) {
    // Definir los rangos de las instancias
    int lower = 1, upper = 30, step = 1;

    // Ejecutar la rutina de pruebas
    std::cerr << "Iniciando verificacion de dataset..." << std::endl << std::endl;

    for (int d = 1; d <= 9; d++) { // Densidades de 0.1 a 0.9
        double density = d / 10.0;

        for (int instance = lower; instance <= upper; instance += step) {
            std::ostringstream oss;
            std::string base_path = "/home/franc/-Sistemas-Adaptativos-MISP/dataset_grafos_no_dirigidos/new_1000_dataset/"; 
            oss << base_path << "erdos_n1000_p0c0.";

            // Manejar los diferentes formatos de densidad
            int density_int = static_cast<int>(density * 100);
            if (density_int % 10 == 0) {
                oss << (density_int / 10);
            } else {
                oss << density_int;
            }

            oss << "_" << instance << ".graph";
            std::string filename = oss.str();
            
            // Llama a la funcion para cargar el grafo y verificar
            // La funcion load_graph_from_file ya tiene la logica para imprimir la verificacion
            load_graph_from_file(filename);

            // Para detener la ejecucion y ver el resultado de una instancia a la vez
            // std::cout << "Presiona ENTER para continuar con el siguiente archivo..." << std::endl;
            // std::cin.ignore();
        }
    }

    std::cerr << std::endl << "Verificacion completada para todos los archivos." << std::endl;
    return 0;
}