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

// Funcion principal que lee el grafo desde un archivo.
void load_graph_from_file(const string &filename) {
    ifstream f(filename);
    if (!f.is_open()) {
        cerr << "Error abriendo archivo: " << filename << "\n";
        n = 0;
        graph.clear();
        return;
    }

    int n_file;
    f >> n_file;
    n = n_file;
    graph.assign(n, Vertice{});

    int u, v;
    int edge_count = 0;
    while (f >> u >> v) {
        if (u >= 0 && u < n && v >= 0 && v < n) {
            graph[u].vecinos.push_back(v);
            graph[v].vecinos.push_back(u);
            edge_count++;
        }
    }
    f.close();

    //  Simple código que verifica que se estén leyendo correctamente los grafos.
    cout << "VERIFICACION DE LECTURA " << endl;
    cout << "Archivo: " << filename << endl;

    // Extraer la densidad del nombre del archivo
    size_t pos_c = filename.find("c0.");
    double file_density = 0.0;
    if (pos_c != string::npos) {
        string density_str = filename.substr(pos_c + 3); // Obtiene la subcadena después de "c"
        size_t pos_underscore = density_str.find('_');
        if (pos_underscore != string::npos) {
            density_str = density_str.substr(0, pos_underscore);
        }
        try {
            file_density = stod("0." + density_str);
        } catch (...) {
            file_density = 0.0; 
        }
    }

    cout << "Numero de vertices: " << n << endl;
    cout << "Numero de aristas: " << edge_count << endl;
    cout << "Densidad del grafo (del nombre del archivo): " << fixed << setprecision(2) << file_density << endl;
    

    cout << "Muestra (primeros 5 vertices):" << endl;
    for(int i = 0; i < min(5, n); ++i) {
        cout << "  Vertice " << i << " tiene " << graph[i].vecinos.size() << " vecinos." << endl;
    }
    cout << "-------------------------------" << endl;
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


