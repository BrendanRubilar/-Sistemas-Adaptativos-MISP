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

struct Vertice {
	std::vector<int> vecinos;
};

std::vector<Vertice> graph;
int n = 0;

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
		std::cerr << "Error al extraer el numero de nodos del nombre del archivo.\n";
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

std::vector<int> local_search(const std::vector<int> &initial_sol) {
	std::vector<int> s = initial_sol;
	bool improved = true;

	while (improved) {
		improved = false;
		std::vector<int> best_neighbor = s;
		int best_size = static_cast<int>(s.size());

		for (int removed_node : s) {
			std::vector<unsigned char> is_in_s(n, 0);
			for (int node : s) {
				is_in_s[node] = 1;
			}
			is_in_s[removed_node] = 0;

			std::vector<int> candidate;
			candidate.reserve(s.size());
			for (int node : s) {
				if (node != removed_node) {
					candidate.push_back(node);
				}
			}

			std::vector<int> potential_adds;
			for (int i = 0; i < n; ++i) {
				if (!is_in_s[i]) {
					bool can_add = true;
					for (int neighbor : graph[i].vecinos) {
						if (is_in_s[neighbor]) {
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
				std::vector<int> temp_sol = candidate;
				temp_sol.push_back(new_node);
				if (static_cast<int>(temp_sol.size()) > best_size) {
					best_neighbor = temp_sol;
					best_size = static_cast<int>(temp_sol.size());
					improved = true;
				}
			}
		}

		if (improved) {
			s = best_neighbor;
		}
	}

	return s;
}

std::vector<int> decode_solution(const std::vector<double> &chromosome) {
	std::vector<int> solution;
	if (chromosome.size() != static_cast<size_t>(n)) {
		return solution;
	}

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

struct Individual {
	std::vector<double> chromosome;
	std::vector<int> solution;
	int fitness = 0;
};

Individual evaluate_chromosome(std::vector<double> chromosome,
							   std::mt19937 &rng,
							   bool use_local_search,
							   double local_search_probability) {
	Individual individual;
	individual.chromosome = std::move(chromosome);
	individual.solution = decode_solution(individual.chromosome);

	if (use_local_search && local_search_probability > 0.0 && !individual.solution.empty()) {
		double draw = std::generate_canonical<double, 10>(rng);
		if (draw <= local_search_probability) {
			std::vector<int> improved = local_search(individual.solution);
			if (improved.size() >= individual.solution.size()) {
				individual.solution = std::move(improved);
			}
		}
	}

	individual.fitness = static_cast<int>(individual.solution.size());
	return individual;
}

int main(int argc, char *argv[]) {
	std::string meta_name = "BRKGA";
	std::string filename = "erdos_n1000_p0c0.1_1.graph";
	int max_seconds = 60;
	int population_size = 60;
	double elite_fraction = 0.25;
	double mutant_fraction = 0.10;
	double inheritance_prob = 0.70;
	bool use_local_search = true;
	double local_search_probability = 0.20;
	unsigned seed = static_cast<unsigned>(std::chrono::high_resolution_clock::now().time_since_epoch().count());

	int arg_index = 1;
	if (arg_index < argc && argv[arg_index][0] != '-') {
		meta_name = argv[arg_index];
		++arg_index;
	}

	for (int i = arg_index; i < argc; ++i) {
		std::string arg = argv[i];
		if (arg == "-i" && i + 1 < argc) {
			filename = argv[++i];
		} else if (arg == "-t" && i + 1 < argc) {
			max_seconds = std::stoi(argv[++i]);
		} else if (arg == "-p" && i + 1 < argc) {
			population_size = std::stoi(argv[++i]);
		} else if (arg == "-e" && i + 1 < argc) {
			elite_fraction = std::stod(argv[++i]);
		} else if (arg == "-m" && i + 1 < argc) {
			mutant_fraction = std::stod(argv[++i]);
		} else if (arg == "-r" && i + 1 < argc) {
			inheritance_prob = std::stod(argv[++i]);
		} else if (arg == "-l" && i + 1 < argc) {
			use_local_search = std::stoi(argv[++i]) != 0;
		} else if (arg == "-L" && i + 1 < argc) {
			local_search_probability = std::stod(argv[++i]);
		} else if (arg == "-s" && i + 1 < argc) {
			seed = static_cast<unsigned>(std::stoul(argv[++i]));
		}
	}

	if (max_seconds <= 0) {
		max_seconds = 1;
	}

	population_size = std::max(4, population_size);
	elite_fraction = std::clamp(elite_fraction, 0.0, 1.0);
	mutant_fraction = std::clamp(mutant_fraction, 0.0, 1.0);
	inheritance_prob = std::clamp(inheritance_prob, 0.0, 1.0);
	local_search_probability = std::clamp(local_search_probability, 0.0, 1.0);

	std::mt19937 rng(seed);

	load_graph(filename);
	if (n == 0) {
		return 1;
	}

	std::cout << "Algoritmo: " << meta_name << " (BRKGA) para MISP\n";
	std::cout << "Instancia: " << filename << "\n";
	std::cout << "Nodos: " << n << "\n";
	std::cout << "Tiempo maximo: " << max_seconds << "s\n";
	std::cout << "Parametros: seed=" << seed
			  << ", poblacion=" << population_size
			  << ", elite=" << elite_fraction
			  << ", mutantes=" << mutant_fraction
			  << ", rho=" << inheritance_prob
			  << ", busqueda_local=" << (use_local_search ? "si" : "no")
			  << ", prob_ls=" << local_search_probability << "\n\n";

	auto start_time = std::chrono::high_resolution_clock::now();

	std::uniform_real_distribution<double> key_dist(0.0, 1.0);

	std::vector<Individual> population;
	population.reserve(population_size);

	std::vector<int> best_solution;
	int best_fitness = 0;
	double best_time = 0.0;

	auto try_update_best = [&](const std::vector<int> &candidate) {
		int candidate_size = static_cast<int>(candidate.size());
		if (candidate_size > best_fitness) {
			auto now = std::chrono::high_resolution_clock::now();
			double elapsed = std::chrono::duration<double>(now - start_time).count();
			best_solution = candidate;
            best_fitness = candidate_size;
            best_time = elapsed;
            std::cout << "Nueva mejor solucion encontrada:" << std::endl;
            std::cout << "Calidad: " << best_fitness << ", Tiempo: " 
                      << std::fixed << std::setprecision(3) << elapsed << "s" << std::endl << std::endl;
		}
	};

	for (int i = 0; i < population_size; ++i) {
		std::vector<double> chromosome(n);
		for (int j = 0; j < n; ++j) {
			chromosome[j] = key_dist(rng);
		}
		Individual individual = evaluate_chromosome(std::move(chromosome), rng, use_local_search, local_search_probability);
		try_update_best(individual.solution);
		population.push_back(std::move(individual));
	}

	bool time_over = false;

	while (!time_over) {
		auto now_loop = std::chrono::high_resolution_clock::now();
		double elapsed_loop = std::chrono::duration<double>(now_loop - start_time).count();
		if (elapsed_loop >= max_seconds) {
			break;
		}

		std::sort(population.begin(), population.end(), [](const Individual &a, const Individual &b) {
			return a.fitness > b.fitness;
		});

		int elite_size = static_cast<int>(std::ceil(population_size * elite_fraction));
		elite_size = std::clamp(elite_size, 1, population_size - 1);

		int mutant_size = static_cast<int>(std::ceil(population_size * mutant_fraction));
		mutant_size = std::max(0, std::min(mutant_size, population_size - elite_size));

		int offspring_count = population_size - elite_size - mutant_size;
		if (offspring_count < 0) {
			mutant_size = std::max(0, mutant_size + offspring_count);
			offspring_count = population_size - elite_size - mutant_size;
		}

		std::vector<Individual> next_population;
		next_population.reserve(population_size);

		for (int i = 0; i < elite_size; ++i) {
			next_population.push_back(population[i]);
		}

		for (int i = 0; i < mutant_size; ++i) {
			auto now_mutant = std::chrono::high_resolution_clock::now();
			double elapsed_mutant = std::chrono::duration<double>(now_mutant - start_time).count();
			if (elapsed_mutant >= max_seconds) {
				time_over = true;
				break;
			}
			std::vector<double> chromosome(n);
			for (int j = 0; j < n; ++j) {
				chromosome[j] = key_dist(rng);
			}
			Individual mutant = evaluate_chromosome(std::move(chromosome), rng, use_local_search, local_search_probability);
			try_update_best(mutant.solution);
			next_population.push_back(std::move(mutant));
		}

		if (time_over) {
			break;
		}

		if (offspring_count > 0) {
			std::uniform_int_distribution<int> elite_dist(0, elite_size - 1);
			std::uniform_int_distribution<int> non_elite_dist(elite_size, population_size - 1);

			for (int i = 0; i < offspring_count; ++i) {
				auto now_child = std::chrono::high_resolution_clock::now();
				double elapsed_child = std::chrono::duration<double>(now_child - start_time).count();
				if (elapsed_child >= max_seconds) {
					time_over = true;
					break;
				}

				const Individual &parent1 = population[elite_dist(rng)];
				const Individual &parent2 = population[non_elite_dist(rng)];

				std::vector<double> chromosome(n);
				for (int j = 0; j < n; ++j) {
					double pick = std::generate_canonical<double, 10>(rng);
					chromosome[j] = (pick < inheritance_prob) ? parent1.chromosome[j] : parent2.chromosome[j];
				}

				Individual child = evaluate_chromosome(std::move(chromosome), rng, use_local_search, local_search_probability);
				try_update_best(child.solution);
				next_population.push_back(std::move(child));
			}
		}

		if (time_over) {
			break;
		}

		if (static_cast<int>(next_population.size()) < population_size) {
			int to_fill = population_size - static_cast<int>(next_population.size());
			for (int i = 0; i < to_fill; ++i) {
				std::vector<double> chromosome(n);
				for (int j = 0; j < n; ++j) {
					chromosome[j] = key_dist(rng);
				}
				Individual filler = evaluate_chromosome(std::move(chromosome), rng, use_local_search, local_search_probability);
				try_update_best(filler.solution);
				next_population.push_back(std::move(filler));
			}
		}

		population = std::move(next_population);
	}
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Fin de la ejecucion. Resultados finales:" << std::endl;
    std::cout << "Mejor solucion encontrada: " << best_fitness << std::endl;
    std::cout << "Tiempo en que se encontro: " << std::fixed << std::setprecision(3) << best_time << "s" << std::endl;
    std::cout << "----------------------------------------" << std::endl;

	return 0;
}
