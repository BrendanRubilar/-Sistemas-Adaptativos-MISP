#include <bits/stdc++.h>
using namespace std;

struct Vertice {
    vector<int> vecinos;
    bool eliminado = false;
};

vector<Vertice> graph;
int n;

void add_edge(int u, int v) {
    graph[u].vecinos.push_back(v);
}

void load_graph_limit_nodes(const string &filename, int max_nodes) {
    ifstream f(filename);
    if (!f.is_open()) {
        cerr << "Error abriendo archivo\n";
        exit(1);
    }

    int n_file;
    f >> n_file;
    
    n = min(max_nodes, n_file);
    graph.assign(n, Vertice{});

    int u, v;
    while (f >> u >> v) {
        if (u < n && v < n) {
            graph[u].vecinos.push_back(v);
            graph[v].vecinos.push_back(u);
        }
    }

    f.close();
}

vector<int> nodos_ordenados_por_grado_simple() {
    vector<pair<int, int>> nodos_con_grado;
    nodos_con_grado.reserve(n);
    
    for (int i = 0; i < n; ++i) {
        int grado = (int)graph[i].vecinos.size();
        nodos_con_grado.emplace_back(grado, i);
    }
    
    // Ordena por grado (ascendente), desempate por ID
    sort(nodos_con_grado.begin(), nodos_con_grado.end());
    
    vector<int> orden;
    orden.reserve(n);
    for (auto& par : nodos_con_grado) {
        orden.push_back(par.second); // agrega el ID del nodo
    }
    
    return orden;
}

vector<int> misp_heuristica_por_orden(const vector<int>& orden) {
    vector<unsigned char> bloqueado(n, 0);
    vector<int> iset;
    iset.reserve(n);

    for (int u : orden) {
        if (!bloqueado[u]) {
            iset.push_back(u);
            bloqueado[u] = 1;
            for (int v : graph[u].vecinos) {
                bloqueado[v] = 1;
            }
        }
    }
    return iset;
}

// int main() {
//     int K = 1000; // cantidad de nodos a cargar (0..K-1)
//     load_graph_limit_nodes("grafo_1.graph", K);

//     // Usa bucket sort para O(n + maxDeg)
//     auto orden = nodos_ordenados_por_grado_simple();
//     auto iset = misp_heuristica_por_orden(orden);

//     cout << "Tamaño IS: " << iset.size() << "\n";
//     for (int u : iset) cout << u << " ";
//     cout << "\n";
//     return 0;
// }
