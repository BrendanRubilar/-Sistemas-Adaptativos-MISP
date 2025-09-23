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
        cerr << "Error abriendo archivo: " << filename << "\n";
        // deja el grafo vacío en caso de fallo
        n = 0;
        graph.clear();
        return;
    }

    string line;
    int n_file = -1;
    vector<pair<int,int>> edges;
    int max_index_in_edges = -1;

    while (getline(f, line)) {
        if (line.empty()) continue;
        // omitir comentarios comunes
        if (line[0] == '#' || line[0] == 'c' || line[0] == '%') continue;

        istringstream ss(line);
        if (n_file == -1) {
            // intentar leer el primer entero como n_file
            if (!(ss >> n_file)) continue;
            // si hay dos enteros en la misma línea, tratar el resto como arista
            int u, v;
            if (ss >> u >> v) {
                edges.emplace_back(u, v);
                max_index_in_edges = max(max_index_in_edges, max(u, v));
            }
        } else {
            int u, v;
            if (ss >> u >> v) {
                edges.emplace_back(u, v);
                max_index_in_edges = max(max_index_in_edges, max(u, v));
            }
        }
    }

    f.close();

    if (n_file <= 0) {
        cerr << "Error: no pude determinar número de nodos en: " << filename << "\n";
        n = 0;
        graph.clear();
        return;
    }

    // Detectar si los índices en el fichero son 1-based
    bool one_based = (max_index_in_edges >= n_file);

    n = min(max_nodes, n_file);
    graph.assign(n, Vertice{});

    for (auto &e : edges) {
        int u = e.first;
        int v = e.second;
        if (one_based) { // convertir a 0-based si es necesario
            --u;
            --v;
        }
        if (u >= 0 && u < n && v >= 0 && v < n) {
            graph[u].vecinos.push_back(v);
            graph[v].vecinos.push_back(u);
        }
        // si la arista tiene nodos fuera del rango solicitado, simplemente se ignora
    }
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

int main(int argc, char *argv[]) {
    string filename = (argc > 1) ? argv[1] : "erdos_n1000_p0c0.1_1.graph";
    int K = (argc > 2) ? stoi(argv[2]) : 1000;

    load_graph_limit_nodes(filename, K);

    cout << "Fichero: " << filename << "\n";
    cout << "Nodos cargados (n): " << n << "\n";

    int sample = min(n, 10);
    cout << "Grados primeros " << sample << " nodos:\n";
    for (int i = 0; i < sample; ++i) {
        cout << i << ": deg=" << graph[i].vecinos.size() << "\n";
    }

    auto orden = nodos_ordenados_por_grado_simple();
    auto iset = misp_heuristica_por_orden(orden);

    cout << "Tamaño IS heurística: " << iset.size() << "\n";
    int show = min((int)iset.size(), 20);
    for (int i = 0; i < show; ++i) {
        cout << iset[i] << (i + 1 == show ? "\n" : " ");
    }

    return 0;
}
