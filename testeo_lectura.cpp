#include <bits/stdc++.h>
#include <filesystem>
using namespace std;
namespace fs = std::filesystem;

struct Vertice {
    vector<int> vecinos;
};

vector<Vertice> graph;
int n;

// -------------------------------
// Función para cargar un grafo desde archivo
// -------------------------------
void load_graph_from_file(const string &filename) {
    ifstream f(filename);
    if (!f.is_open()) {
        cerr << "Error abriendo archivo: " << filename << "\n";
        exit(1);
    }

    int n_file;
    f >> n_file; // número de nodos indicado en el archivo
    n = n_file;
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

// -------------------------------
// Extraer densidad desde el nombre del archivo
// -------------------------------
double extract_density_from_filename(const string &filename) {
    size_t cpos = filename.find("c");
    if (cpos == string::npos) return -1.0;

    size_t underscore = filename.find("_", cpos);
    string densidad_str = filename.substr(cpos + 1, underscore - cpos - 1);
    return stod(densidad_str);
}

// -------------------------------
// MAIN
// -------------------------------
int main() {
    string base_path = "dataset_grafos_no_dirigidos/new_1000_dataset";

    // 🔹 ESTE for debe estar dentro del main
    for (auto &entry : fs::directory_iterator(base_path)) {
        string filepath = entry.path().string();
        string filename = entry.path().filename().string();

        // Extraer densidad
        double densidad = extract_density_from_filename(filename);

        // Cargar el grafo
        load_graph_from_file(filepath);

        // Mostrar información básica
        cout << "Archivo: " << filename << "\n";
        cout << " - Nodos: " << n << "\n";
        cout << " - Densidad (extraída del nombre): " << densidad << "\n";

        long long total_edges = 0;
        for (int i = 0; i < n; i++) total_edges += graph[i].vecinos.size();
        cout << " - Aristas cargadas: " << total_edges / 2 << "\n\n";
    }

    return 0;
}
