#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>

// Forward declarations for the shared functions
std::vector<int> nodos_ordenados_por_grado_simple();
std::vector<int> misp_heuristica_por_orden(const std::vector<int>& orden);
void load_graph_limit_nodes(const std::string& filename, int max_nodes);

struct Vertice {
    std::vector<int> vecinos;
    bool eliminado = false;
};

extern std::vector<Vertice> graph;
extern int n;

#endif // GRAPH_H