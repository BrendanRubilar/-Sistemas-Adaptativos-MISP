#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <random>
#include <chrono>
#include <iomanip>
#include <regex>

// Estructura de Vertice y variables globales
struct Vertice {
    std::vector<int> vecinos;
};

std::vector<Vertice> graph;
int n;

// Funciones para cargar el grafo y otras utilidades (misma implementación anterior)
void load_graph(const std::string &filename) {
    std::ifstream f(filename);
    if (!f.is_open()) {
        std::cerr << "Error abriendo archivo: " << filename << "\n";
        n = 0;
        graph.clear();
        return;
    }

    std::string line;
    int n_file = -1;
    std::vector<std::pair<int,int>> edges;
    int max_index_in_edges = -1;

    // Se extrae el número de nodos del nombre del archivo
    std::regex node_regex("n([0-9]+)");
    std::smatch matches;
    if (std::regex_search(filename, matches, node_regex) && matches.size() > 1) {
        n_file = std::stoi(matches[1].str());
    } else {
        std::cerr << "Error: No se pudo extraer el numero de nodos del nombre del archivo.\n";
        n = 0;
        graph.clear();
        return;
    }

    while (getline(f, line)) {
        if (line.empty() || line[0] == '#' || line[0] == 'c' || line[0] == '%') continue;
        std::istringstream ss(line);
        int u, v;
        if (ss >> u >> v) {
            edges.emplace_back(u, v);
            max_index_in_edges = std::max(max_index_in_edges, std::max(u, v));
        }
    }

    f.close();

    bool one_based = (max_index_in_edges >= n_file);
    n = n_file;
    graph.assign(n, Vertice{});

    for (auto &e : edges) {
        int u = e.first;
        int v = e.second;
        if (one_based) {
            --u;
            --v;
        }
        if (u >= 0 && u < n && v >= 0 && v < n) {
            graph[u].vecinos.push_back(v);
            graph[v].vecinos.push_back(u);
        }
    }
}

// Fase de Búsqueda Local de GRASP (misma implementación anterior)
std::vector<int> local_search(const std::vector<int>& initial_solution) {
    std::vector<int> s = initial_solution;
    bool improved = true;

    while (improved) {
        improved = false;
        std::vector<int> best_neighbor_solution = s;
        int best_size = s.size();

        for (int removed_node : s) {
            std::vector<unsigned char> is_in_s(n, 0);
            for(int node : s) is_in_s[node] = 1;
            is_in_s[removed_node] = 0;

            std::vector<int> current_neighbor_solution;
            for(int node : s) {
                if (node != removed_node) {
                    current_neighbor_solution.push_back(node);
                }
            }

            std::vector<int> potential_additions;
            for(int i = 0; i < n; ++i) {
                if (!is_in_s[i]) {
                    bool can_add = true;
                    for (int neighbor : graph[i].vecinos) {
                        if (is_in_s[neighbor]) {
                            can_add = false;
                            break;
                        }
                    }
                    if (can_add) {
                        potential_additions.push_back(i);
                    }
                }
            }
            
            for (int new_node : potential_additions) {
                std::vector<int> temp_solution = current_neighbor_solution;
                temp_solution.push_back(new_node);

                if (temp_solution.size() > best_size) {
                    best_neighbor_solution = temp_solution;
                    best_size = temp_solution.size();
                    improved = true;
                }
            }
        }
        if (improved) {
            s = best_neighbor_solution;
        }
    }
    return s;
}

// Fase de Construcción Aleatorizada de GRASP (misma implementación anterior)
std::vector<int> build_randomized_greedy_solution(double alpha, std::mt19937& g) {
    std::vector<unsigned char> is_available(n, 1);
    std::vector<int> iset;
    std::vector<int> available_nodes;
    for(int i = 0; i < n; ++i) {
        available_nodes.push_back(i);
    }

    while (!available_nodes.empty()) {
        std::vector<std::pair<int, int>> candidates_with_degree;
        int max_degree = 0;
        int min_degree = n;
        
        for (int node : available_nodes) {
            int degree = 0;
            for (int neighbor : graph[node].vecinos) {
                if (is_available[neighbor]) {
                    degree++;
                }
            }
            candidates_with_degree.push_back({degree, node});
            max_degree = std::max(max_degree, degree);
            min_degree = std::min(min_degree, degree);
        }

        std::vector<int> rcl;
        double threshold = min_degree + alpha * (max_degree - min_degree);
        
        for (const auto& p : candidates_with_degree) {
            if (p.first <= threshold) {
                rcl.push_back(p.second);
            }
        }

        if (rcl.empty()) break;
        
        std::uniform_int_distribution<> dis(0, rcl.size() - 1);
        int u = rcl[dis(g)];
        
        iset.push_back(u);
        is_available[u] = 0;
        
        for (int v : graph[u].vecinos) {
            is_available[v] = 0;
        }

        available_nodes.clear();
        for(int i = 0; i < n; ++i) {
            if(is_available[i]) {
                available_nodes.push_back(i);
            }
        }
    }
    return iset;
}

std::vector<int> run_grasp(const std::string& filename, int max_seconds, double alpha, std::uint64_t seed) {
    // Configuración de la semilla para la aleatoriedad
    std::mt19937 g(seed);

    // Cargar el grafo
    load_graph(filename);

    if (n == 0) {
        // Retorna un vector vacío para indicar un error en la carga del grafo
        return {}; 
    }

    std::cout << "Algoritmo: GRASP para MISP\n";
    std::cout << "Instancia: " << filename << "\n";
    std::cout << "Nodos: " << n << "\n";
    std::cout << "Tiempo maximo: " << max_seconds << "s\n";
    std::cout << "Parametros: alpha=" << alpha << ", seed=" << seed << "\n\n";

    std::vector<int> best_solution;
    double best_time = 0.0;
    
    auto start_time = std::chrono::high_resolution_clock::now();

    while (true) {
        auto current_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = current_time - start_time;
        if (elapsed.count() >= max_seconds) {
            break;
        }

        std::vector<int> current_solution = build_randomized_greedy_solution(alpha, g);
        std::vector<int> improved_solution = local_search(current_solution);
        
        if (improved_solution.size() > best_solution.size()) {
            best_solution = improved_solution;
            best_time = elapsed.count();
            // Salida para el Any-Time Behavior
            std::cout << "Nueva mejor solucion encontrada:\n";
            std::cout << "Calidad: " << best_solution.size() << ", Tiempo: " << std::fixed << std::setprecision(3) << best_time << "s\n\n";
        }
    }

    // Salida final
    std::cout << "----------------------------------------\n";
    std::cout << "Fin de la ejecucion. Resultados finales:\n";
    std::cout << "Mejor solucion encontrada: " << best_solution.size() << "\n";
    std::cout << "Tiempo en que se encontro: " << std::fixed << std::setprecision(3) << best_time << "s\n";
    std::cout << "----------------------------------------\n";

    return best_solution;
}

// Función principal adaptada para cumplir con los requisitos
/*int main(int argc, char *argv[]) {
    // Parámetros por defecto
    std::string filename = "erdos_n1000_p0c0.1_1.graph";
    int max_seconds = 60;
    double alpha = 0.5;
    unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();

    // Análisis de los argumentos de la línea de comandos
    if (argc > 1) {
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "-i" && i + 1 < argc) {
                filename = argv[++i];
            } else if (arg == "-t" && i + 1 < argc) {
                max_seconds = std::stoi(argv[++i]);
            } else if (arg == "-a" && i + 1 < argc) {
                alpha = std::stod(argv[++i]);
            } else if (arg == "-s" && i + 1 < argc) {
                seed = std::stoul(argv[++i]);
            }
        }
    }
    
    // Configuración de la semilla para la aleatoriedad
    std::mt19937 g(seed);

    // Cargar el grafo sin un límite explícito de nodos, se obtiene del nombre del archivo
    load_graph(filename);

    if (n == 0) {
        return 1;
    }

    std::cout << "Algoritmo: GRASP para MISP\n";
    std::cout << "Instancia: " << filename << "\n";
    std::cout << "Nodos: " << n << "\n";
    std::cout << "Tiempo maximo: " << max_seconds << "s\n";
    std::cout << "Parametros: alpha=" << alpha << ", seed=" << seed << "\n\n";

    std::vector<int> best_solution;
    double best_time = 0.0;
    
    auto start_time = std::chrono::high_resolution_clock::now();

    while (true) {
        auto current_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = current_time - start_time;
        if (elapsed.count() >= max_seconds) {
            break;
        }

        std::vector<int> current_solution = build_randomized_greedy_solution(alpha, g);
        std::vector<int> improved_solution = local_search(current_solution);
        
        if (improved_solution.size() > best_solution.size()) {
            best_solution = improved_solution;
            best_time = elapsed.count();
            // Salida para el Any-Time Behavior
            std::cout << "Nueva mejor solucion encontrada:\n";
            std::cout << "Calidad: " << best_solution.size() << ", Tiempo: " << std::fixed << std::setprecision(3) << best_time << "s\n\n";
        }
    }

    // Salida final
    std::cout << "----------------------------------------\n";
    std::cout << "Fin de la ejecucion. Resultados finales:\n";
    std::cout << "Mejor solucion encontrada: " << best_solution.size() << "\n";
    std::cout << "Tiempo en que se encontro: " << std::fixed << std::setprecision(3) << best_time << "s\n";
    std::cout << "----------------------------------------\n";

    return 0;
}
*/