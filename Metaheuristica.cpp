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


struct Vertice {
    std::vector<int> vecinos;
};

std::vector<Vertice> graph;
int n;

void load_graph(const std::string &filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error abriendo archivo: " << filename << "\n";
        n = 0;
        graph.clear();
        return;
    }

    std::string line;
    int n_from_file = -1;
    std::vector<std::pair<int,int>> edges;
    int max_idx = -1;

    std::regex node_regex("n([0-9]+)");
    std::smatch matches;
    if (std::regex_search(filename, matches, node_regex) && matches.size() > 1) {
        n_from_file = std::stoi(matches[1].str());
    } else {
        std::cerr << "Error al extraer el número de nodos del nombre del archivo.\n";
        n = 0;
        graph.clear();
        return;
    }

    while (getline(file, line)) {
        if (line.empty() || line[0] == '#' || line[0] == 'c' || line[0] == '%') continue;
        std::istringstream ss(line);
        int u, v;
        if (ss >> u >> v) {
            edges.emplace_back(u, v);
            max_idx = std::max(max_idx, std::max(u, v));
        }
    }

    file.close();

    bool one_based = (max_idx >= n_from_file);
    n = n_from_file;
    graph.assign(n, Vertice{});

    for (auto &e : edges) {
        int u = e.first;
        int v = e.second;
        if (one_based) {
            --u;
            --v;
        }
        if (u>=0 && u<n && v>=0 && v<n) {
            graph[u].vecinos.push_back(v);
            graph[v].vecinos.push_back(u);
        }
    }
}

std::vector<int> local_search(const std::vector<int>& initial_sol) {
    std::vector<int> s = initial_sol;
    bool improved = true;

    while (improved) {
        improved = false;
        std::vector<int> best_vecino_sol = s;
        int best_size = s.size();

        for (int removed_node : s) {
            std::vector<unsigned char> is_in_s(n, 0);
            for(int node : s) is_in_s[node] = 1;
            is_in_s[removed_node] = 0;

            std::vector<int> actual_vecino_sol;
            for(int node : s) {
                if (node != removed_node) {
                    actual_vecino_sol.push_back(node);
                }
            }

            std::vector<int> potential_adds;
            for(int i = 0; i < n; ++i) {
                if (!is_in_s[i]) {
                    bool can_add = true;
                    for (int vecino : graph[i].vecinos) {
                        if (is_in_s[vecino]) {
                            can_add = false;
                            break;
                        }
                    }
                    if (can_add) {
                        potential_adds.push_back(i);
                    }
                }
            }
            
            for (int new_node : potential_adds) {
                std::vector<int> temp_sol = actual_vecino_sol;
                temp_sol.push_back(new_node);

                if (temp_sol.size() > best_size) {
                    best_vecino_sol = temp_sol;
                    best_size = temp_sol.size();
                    improved = true;
                }
            }
        }
        if (improved) {
            s = best_vecino_sol;
        }
    }
    return s;
}

std::vector<int> build_randomized_greedy_solution(double alpha, std::mt19937& g) {
    std::vector<unsigned char> available_status(n, 1);
    std::vector<int> iset;
    std::vector<int> avail_nodes;
    for(int i = 0; i < n; ++i) {
        avail_nodes.push_back(i);
    }

    while (!avail_nodes.empty()) {
        std::vector<std::pair<int, int>> candidates_deg;
        int max_deg = 0;
        int min_deg = n;
        
        for (int node : avail_nodes) {
            int degree = 0;
            for (int vecino : graph[node].vecinos) {
                if (available_status[vecino]) {
                    degree++;
                }
            }
            candidates_deg.push_back({degree, node});
            max_deg = std::max(max_deg, degree);
            min_deg = std::min(min_deg, degree);
        }

        std::vector<int> rcl;
        double limit = min_deg + alpha * (max_deg - min_deg);
        
        for (const auto& p : candidates_deg) {
            if (p.first <= limit) {
                rcl.push_back(p.second);
            }
        }

        if (rcl.empty()) break;
        
        std::uniform_int_distribution<> dis(0, rcl.size() - 1);
        int u = rcl[dis(g)];
        
        iset.push_back(u);
        available_status[u] = 0;
        
        for (int v : graph[u].vecinos) {
            available_status[v] = 0;
        }

        avail_nodes.clear();
        for(int i = 0; i < n; ++i) {
            if(available_status[i]) {
                avail_nodes.push_back(i);
            }
        }
    }
    return iset;
}

std::vector<int> run_grasp(const std::string& filename, int max_seconds, double alpha, std::uint32_t seed) {
    std::mt19937 g(seed);

    load_graph(filename);

    if (n == 0) { // En caso de error, no retornar nada
        return {}; 
    }

    std::cout << "Algoritmo: GRASP para MISP\n";
    std::cout << "Instancia: " << filename << "\n";
    std::cout << "Nodos: " << n << "\n";
    std::cout << "Tiempo maximo: " << max_seconds << "s\n";
    std::cout << "Parametros: alpha=" << alpha << ", seed=" << seed << "\n\n";

    std::vector<int> best_sol;
    double best_time = 0.0;
    
    auto start_time = std::chrono::high_resolution_clock::now();

    while (true) {
        auto current_time = std::chrono::high_resolution_clock::now(); 
        std::chrono::duration<double> elapsed = current_time - start_time; 
        if (elapsed.count() >= max_seconds) {
            break;
        }

        std::vector<int> current_sol = build_randomized_greedy_solution(alpha, g);
        std::vector<int> improved_sol = local_search(current_sol);
        
        if (improved_sol.size() > best_sol.size()) {
            best_sol = improved_sol;
            best_time = elapsed.count();

            std::cout << "Nueva mejor solucion encontrada:\n";
            std::cout << "Calidad: " << best_sol.size() << ", Tiempo: " << std::fixed << std::setprecision(3) << best_time << "s\n\n";
        }
    }

    // Salida final
    std::cout << "----------------------------------------\n";
    std::cout << "Fin de la ejecucion. Resultados finales:\n";
    std::cout << "Mejor solucion encontrada: " << best_sol.size() << "\n";
    std::cout << "Tiempo en que se encontro: " << std::fixed << std::setprecision(3) << best_time << "s\n";
    std::cout << "----------------------------------------\n";

    return best_sol;
}

// Función principal adaptada para cumplir con los requisitos
int main(int argc, char *argv[]) {
    // Parámetros por defecto
    std::string filename = "erdos_n1000_p0c0.1_1.graph";
    int max_seconds = 60;
    double alpha = 0.5;
    unsigned seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();

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

    std::mt19937 g(seed);

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

        std::vector<int> actual_solution = build_randomized_greedy_solution(alpha, g);
        std::vector<int> improved_solution = local_search(actual_solution);
        
        if (improved_solution.size() > best_solution.size()) {
            best_solution = improved_solution;
            best_time = elapsed.count();
            std::cout << "Nueva mejor solucion encontrada:\n";
            std::cout << "Calidad: " << best_solution.size() << ", Tiempo: " << std::fixed << std::setprecision(3) << best_time << "s\n\n";
        }
    }

    std::cout << "----------------------------------------\n";
    std::cout << "Fin de la ejecucion. Resultados finales:\n";
    std::cout << "Mejor solucion encontrada: " << best_solution.size() << "\n";
    std::cout << "Tiempo en que se encontro: " << std::fixed << std::setprecision(3) << best_time << "s\n";
    std::cout << "----------------------------------------\n";

    return 0;
}
