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
#include <numeric>
#include <cmath>

// Definición de tipos de tiempo para compatibilidad con C++17 sin extensions
using TimePoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

// --- I. Definiciones Globales y Estructuras ---
struct Vertice {
    std::vector<int> vecinos;
};

std::vector<Vertice> graph;
int n = 0;

// Estructura para encapsular todos los hiperparámetros
struct Params {
    std::string meta_name;
    std::string filename;
    int max_seconds;
    int population_size;
    double elite_fraction;
    double mutant_fraction;
    double inheritance_prob;
    double alpha_grasp; 
    bool use_local_search;
    double local_search_probability;
    unsigned seed;
};

struct Individual {
    std::vector<double> chromosome;
    std::vector<int> solution;
    int fitness = 0;
};

// --- DECLARACIONES DE FUNCIONES (Para modularidad) ---
void try_update_best(const std::vector<int> &candidate, const TimePoint& start_time,
                     double& best_time, int& best_fitness, std::vector<int>& best_solution);
std::vector<int> local_search(const std::vector<int> &initial_sol);
std::vector<int> build_randomized_greedy_solution(double alpha, std::mt19937& g);
Individual evaluate_chromosome(std::vector<double> chromosome,
                              std::mt19937 &rng,
                              const Params& params,
                              bool is_mutant);
// --- FIN DECLARACIONES ---


// --- II. Funciones de Carga y Utilidades ---

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
    std::vector<std::pair<int, int>> edges;
    int max_idx = -1;

    std::regex node_regex("n([0-9]+)");
    std::smatch matches;
    if (std::regex_search(filename, matches, node_regex) && matches.size() > 1) {
        n_from_file = std::stoi(matches[1].str());
    } else {
        // Asumimos 0 si no se puede leer el N de nodos
        n_from_file = 0; 
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

    for (auto &edge : edges) {
        int u = edge.first;
        int v = edge.second;
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


// Función de utilidad para actualizar la mejor solución (CORREGIDA LA FIRMA)
void try_update_best(const std::vector<int> &candidate, const TimePoint& start_time,
                     double& best_time, int& best_fitness, std::vector<int>& best_solution) {
    int candidate_size = static_cast<int>(candidate.size());
    if (candidate_size > best_fitness) {
        auto now = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(now - start_time).count();
        best_solution = candidate;
        best_fitness = candidate_size;
        best_time = elapsed;
        std::cout << "Nueva mejor solucion encontrada:\n";
        std::cout << "Calidad: " << best_fitness << ", Tiempo: " 
                  << std::fixed << std::setprecision(3) << elapsed << "s\n\n";
    }
}


// --- III. Núcleo GRASP y Búsqueda Local ---

std::vector<int> local_search(const std::vector<int> &initial_sol) {
    // [Implementación de la Búsqueda Local (LS)]
    std::vector<int> s = initial_sol;
    bool improved = true;

    while (improved) {
        improved = false;
        std::vector<int> best_neighbor = s;
        int best_size = static_cast<int>(s.size());

        for (int removed_node : s) {
            std::vector<unsigned char> is_in_s(n, 0);
            for (int node : s) is_in_s[node] = 1;
            is_in_s[removed_node] = 0;

            std::vector<int> candidate;
            candidate.reserve(s.size());
            for (int node : s) {
                if (node != removed_node) candidate.push_back(node);
            }

            std::vector<int> potential_adds;
            for (int i = 0; i < n; ++i) {
                if (!is_in_s[i]) {
                    bool can_add = true;
                    for (int neighbor : graph[i].vecinos) {
                        if (is_in_s[neighbor]) { can_add = false; break; }
                    }
                    if (can_add) potential_adds.push_back(i);
                }
            }

            for (int new_node : potential_adds) {
                std::vector<int> temp_sol = candidate;
                temp_sol.push_back(new_node);
                if (static_cast<int>(temp_sol.size()) > best_size) {
                    best_neighbor = temp_sol;
                    best_size = static_cast<int>(temp_sol.size());
                    improved = true;
                    // Primera mejora (First-Improvement) para acelerar
                    s = best_neighbor;
                    goto end_local_search_iter; 
                }
            }
        }
        end_local_search_iter:;

        if (improved) {
            s = best_neighbor;
        }
    }
    return s;
}

std::vector<int> build_randomized_greedy_solution(double alpha, std::mt19937& g) {
    // [Implementación de la Construcción Aleatoria (GRASP)]
    std::vector<unsigned char> available_status(n, 1);
    std::vector<int> iset;
    std::vector<int> avail_nodes;
    for(int i = 0; i < n; ++i) avail_nodes.push_back(i);

    while (!avail_nodes.empty()) {
        std::vector<std::pair<int, int>> candidates_deg;
        int max_deg = 0;
        int min_deg = n;
        
        for (int node : avail_nodes) {
            int degree = 0;
            for (int vecino : graph[node].vecinos) {
                if (available_status[vecino]) degree++;
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

// --- IV. Núcleo BRKGA ---

std::vector<int> decode_solution(const std::vector<double> &chromosome) {
    // [Decodificador BRKGA]
    std::vector<int> solution;
    if (chromosome.size() != static_cast<size_t>(n)) return solution;

    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return chromosome[a] < chromosome[b];
    });

    std::vector<unsigned char> blocked(n, 0);
    for (int node : order) {
        if (!blocked[node]) {
            solution.push_back(node);
            blocked[node] = 1;
            for (int neighbor : graph[node].vecinos) {
                blocked[neighbor] = 1;
            }
        }
    }
    return solution;
}

Individual evaluate_chromosome(std::vector<double> chromosome,
                              std::mt19937 &rng,
                              const Params& params,
                              bool is_mutant) {
    
    Individual individual;
    
    if (is_mutant) {
        // Hibridación: Los mutantes se generan con la Construcción Aleatoria GRASP
        individual.solution = build_randomized_greedy_solution(params.alpha_grasp, rng);
        individual.chromosome.assign(n, 0.0);
        
    } else {
        // Individuo decodificado del BRKGA (hijo de crossover o inicial)
        individual.chromosome = std::move(chromosome);
        individual.solution = decode_solution(individual.chromosome);
    }

    // Refinamiento Vigoroso (Búsqueda Local) - Aplicado a Mutantes y Descendientes BRKGA
    if (params.use_local_search && !individual.solution.empty()) {
        double draw = std::generate_canonical<double, 10>(rng);
        if (draw <= params.local_search_probability) {
            std::vector<int> improved = local_search(individual.solution);
            if (improved.size() > individual.solution.size()) {
                individual.solution = std::move(improved);
            }
        }
    }

    individual.fitness = static_cast<int>(individual.solution.size());
    return individual;
}

// --- V. Funciones Modulares de Control ---

Params parse_arguments(int argc, char *argv[]) {
    // Parámetros por Defecto (Óptimo BRKGA)
    Params p;
    p.meta_name = "HIBRIDO-BRKGA";
    p.filename = "erdos_n1000_p0c0.1_1.graph";
    p.max_seconds = 60;
    p.population_size = 35;       
    p.elite_fraction = 0.39;   
    p.mutant_fraction = 0.15;  
    p.inheritance_prob = 0.61; 
    p.use_local_search = true;
    p.local_search_probability = 0.20; 
    p.alpha_grasp = 0.20;              
    p.seed = static_cast<unsigned>(std::chrono::high_resolution_clock::now().time_since_epoch().count());

    int arg_index = 1;
    if (arg_index < argc && argv[arg_index][0] != '-') {
        p.meta_name = argv[arg_index];
        ++arg_index;
    }

    for (int i = arg_index; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-i" && i + 1 < argc) { p.filename = argv[++i]; } 
        else if (arg == "-t" && i + 1 < argc) { p.max_seconds = std::stoi(argv[++i]); } 
        else if (arg == "-p" && i + 1 < argc) { p.population_size = std::stoi(argv[++i]); } 
        else if (arg == "-e" && i + 1 < argc) { p.elite_fraction = std::stod(argv[++i]); } 
        else if (arg == "-m" && i + 1 < argc) { p.mutant_fraction = std::stod(argv[++i]); } 
        else if (arg == "-r" && i + 1 < argc) { p.inheritance_prob = std::stod(argv[++i]); } 
        else if (arg == "-l" && i + 1 < argc) { p.use_local_search = std::stoi(argv[++i]) != 0; } 
        else if (arg == "-L" && i + 1 < argc) { p.local_search_probability = std::stod(argv[++i]); } 
        else if (arg == "-a" && i + 1 < argc) { p.alpha_grasp = std::stod(argv[++i]); } 
        else if (arg == "-s" && i + 1 < argc) { p.seed = static_cast<unsigned>(std::stoul(argv[++i])); }
    }

    // Validación de parámetros (uso de min/max como alternativa robusta a std::clamp)
    if (p.max_seconds <= 0) p.max_seconds = 1;
    p.population_size = std::max(4, p.population_size);
    p.elite_fraction = std::max(0.0, std::min(p.elite_fraction, 1.0));
    p.mutant_fraction = std::max(0.0, std::min(p.mutant_fraction, 1.0));
    p.inheritance_prob = std::max(0.0, std::min(p.inheritance_prob, 1.0));
    p.local_search_probability = std::max(0.0, std::min(p.local_search_probability, 1.0));
    p.alpha_grasp = std::max(0.0, std::min(p.alpha_grasp, 1.0));
    
    return p;
}

void initialize_population(std::vector<Individual>& population, const Params& params, 
                         std::mt19937& rng, double& best_time, int& best_fitness, 
                         std::vector<int>& best_solution, const TimePoint& start_time) {
    
    std::uniform_real_distribution<double> key_dist(0.0, 1.0);
    population.reserve(params.population_size);

    for (int i = 0; i < params.population_size; ++i) {
        std::vector<double> chromosome(n);
        for (int j = 0; j < n; ++j) {
            chromosome[j] = key_dist(rng);
        }
        // Inicialización: No son mutantes
        Individual individual = evaluate_chromosome(std::move(chromosome), rng, params, false);
        try_update_best(individual.solution, start_time, best_time, best_fitness, best_solution);
        population.push_back(std::move(individual));
    }
}

void run_brkga_generation(std::vector<Individual>& population, const Params& params, 
                         std::mt19937& rng, const TimePoint& start_time, 
                         double& best_time, int& best_fitness, std::vector<int>& best_solution,
                         bool& time_over) {
    
    // Si ya se acabó el tiempo, salir.
    if (std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count() >= params.max_seconds) {
        time_over = true;
        return;
    }
    
    std::uniform_real_distribution<double> key_dist(0.0, 1.0);

    // 1. Clasificación y Partición de la Población
    std::sort(population.begin(), population.end(), [](const Individual &a, const Individual &b) {
        return a.fitness > b.fitness;
    });

    int elite_size = static_cast<int>(std::ceil(params.population_size * params.elite_fraction));
    elite_size = std::max(1, std::min(elite_size, params.population_size - 1));

    int mutant_size = static_cast<int>(std::ceil(params.population_size * params.mutant_fraction));
    mutant_size = std::max(0, std::min(mutant_size, params.population_size - elite_size));

    int offspring_count = params.population_size - elite_size - mutant_size;
    
    std::vector<Individual> next_population;
    next_population.reserve(params.population_size);

    // A. Elitismo
    for (int i = 0; i < elite_size; ++i) {
        next_population.push_back(population[i]);
    }

    // B. Mutantes Híbridos (Construcción GRASP)
    for (int i = 0; i < mutant_size; ++i) {
        if (std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count() >= params.max_seconds) { time_over = true; break; }
        
        // Llama a evaluate_chromosome con is_mutant=true para activar la construcción GRASP
        Individual mutant = evaluate_chromosome(std::vector<double>(n), rng, params, true); 
        try_update_best(mutant.solution, start_time, best_time, best_fitness, best_solution);
        next_population.push_back(std::move(mutant));
    }
    if (time_over) return;

    // C. Descendientes por Cruce (Crossover BRKGA)
    if (offspring_count > 0) {
        std::uniform_int_distribution<int> elite_dist(0, elite_size - 1);
        std::uniform_int_distribution<int> non_elite_dist(elite_size, params.population_size - 1);

        for (int i = 0; i < offspring_count; ++i) {
            if (std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_time).count() >= params.max_seconds) { time_over = true; break; }

            const Individual &parent1 = population[elite_dist(rng)];
            const Individual &parent2 = population[non_elite_dist(rng)];

            // Crossover Sesgado
            std::vector<double> chromosome(n);
            for (int j = 0; j < n; ++j) {
                double pick = std::generate_canonical<double, 10>(rng);
                chromosome[j] = (pick < params.inheritance_prob) ? parent1.chromosome[j] : parent2.chromosome[j];
            }

            // Evalúa y aplica Búsqueda Local (LS) sobre el descendiente decodificado
            Individual child = evaluate_chromosome(std::move(chromosome), rng, params, false);
            try_update_best(child.solution, start_time, best_time, best_fitness, best_solution);
            next_population.push_back(std::move(child));
        }
    }
    if (time_over) return;

    // 2. Reemplazo y Actualización
    if (static_cast<int>(next_population.size()) == params.population_size) {
         population = std::move(next_population);
    }
}


// --- VI. MAIN Modular ---
int main(int argc, char *argv[]) {
    // 1. Lectura de Parámetros
    Params params = parse_arguments(argc, argv);
    std::mt19937 rng(params.seed);

    // 2. Carga de Grafo
    load_graph(params.filename);
    if (n == 0) return 1;

    // 3. Variables de Seguimiento (Best Solution)
    std::vector<Individual> population;
    std::vector<int> best_solution;
    int best_fitness = 0;
    double best_time = 0.0;
    
    // CORREGIDO: Usamos TimePoint en lugar de auto
    TimePoint start_time = std::chrono::high_resolution_clock::now();

    // Impresión de parámetros iniciales
    std::cout << "Algoritmo: " << params.meta_name << " (Híbrido BRKGA-GRASP) para MISP\n";
    std::cout << "Instancia: " << params.filename << ", Nodos: " << n << "\n";
    std::cout << "Tiempo maximo: " << params.max_seconds << "s\n";
    std::cout << "Parametros: seed=" << params.seed
              << ", P=" << params.population_size
              << ", e=" << params.elite_fraction
              << ", m_GRASP=" << params.mutant_fraction
              << ", rho=" << params.inheritance_prob
              << ", alpha=" << params.alpha_grasp 
              << ", prob_ls=" << params.local_search_probability << "\n\n";

    // 4. Inicialización
    initialize_population(population, params, rng, best_time, best_fitness, best_solution, start_time);

    // 5. Ciclo Principal (Modular)
    bool time_over = false;
    while (!time_over) {
        // Ejecutar una generación modular
        run_brkga_generation(population, params, rng, start_time, best_time, best_fitness, best_solution, time_over);
    }
    
    // 6. Salida Final
    TimePoint end_time = std::chrono::high_resolution_clock::now();
    double total_elapsed = std::chrono::duration<double>(end_time - start_time).count();
    
    std::cout << "\n----------------------------------------\n";
    std::cout << "Fin de la ejecucion. Resultados finales:\n";
    std::cout << "Mejor solucion encontrada: " << best_fitness << "\n";
    std::cout << "Tiempo en que se encontro: " << std::fixed << std::setprecision(3) << best_time << "s\n";
    std::cout << "Tiempo total de ejecucion: " << std::fixed << std::setprecision(3) << total_elapsed << "s\n";
    std::cout << "----------------------------------------\n";
    
    // Salida final para irace (Coste Negativo)
    std::cout << "-" << best_fitness << " " << total_elapsed << std::endl;

    return 0;
}