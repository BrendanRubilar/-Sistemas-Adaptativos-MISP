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

void load_graph(const string &filename) {
    ifstream f(filename);
    if (!f.is_open()) {
        cerr << "Error abriendo archivo: " << filename << "\n";
        n = 0;
        graph.clear();
        return;
    }

    f >> n;
    graph.assign(n, Vertice{});

    int u, v;
    while (f >> u >> v) {
        if (u >= 0 && u < n && v >= 0 && v < n) {
            graph[u].vecinos.push_back(v);
            graph[v].vecinos.push_back(u); // grafo no dirigido
        }
    }

    f.close();
}

vector<int> misp_heuristica_greedy() {
    vector<unsigned char> bloqueado(n, 0);
    vector<int> iset;
    iset.reserve(n);

    vector<pair<int, int>> nodos_con_grado;
    nodos_con_grado.reserve(n);
    
    for (int i = 0; i < n; ++i) {
        int grado = (int)graph[i].vecinos.size();
        nodos_con_grado.emplace_back(grado, i);
    }
    
    sort(nodos_con_grado.begin(), nodos_con_grado.end());

    for (auto& par : nodos_con_grado) {
        int u = par.second;
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

vector<int> misp_greedy_aleatorizado(float alpha) {
    vector<unsigned char> disponible(n, 1);
    vector<int> iset;
    iset.reserve(n);

    int restantes = n;

    while (restantes > 0) {
        vector<pair<int,int>> candidatos;
        candidatos.reserve(restantes);
        for (int u = 0; u < n; ++u) {
            if (disponible[u]) {
                int grado = (int)graph[u].vecinos.size();
                candidatos.emplace_back(grado, u);
            }
        }

        if (candidatos.empty()) break;

        sort(candidatos.begin(), candidatos.end());

        int tam_bloque = max(1, (int)ceil(candidatos.size() * alpha));

        int idx = rand() % tam_bloque;
        int elegido = candidatos[idx].second;

        iset.push_back(elegido);

        disponible[elegido] = 0;
        restantes--;
        for (int v : graph[elegido].vecinos) {
            if (disponible[v]) {
                disponible[v] = 0;
                restantes--;
            }
        }
    }

    return iset;
}

int main2(int argc, char* argv[]) {
    srand(time(nullptr));
    auto print_usage = [](){
        cerr << "Uso:\n"
             << "  ./Heuristica Greedy -i <archivo_instancia>\n"
             << "  ./Heuristica Greedy-probabilista -i <archivo_instancia> <alpha>\n"
             << "    alpha en (0,1], fracción de candidatos considerados.\n";
    };

    if (argc < 4) {
        print_usage();
        return 1;
    }

    string metodo = argv[1];
    string flag = argv[2];
    if (flag != "-i") {
        cerr << "Error: falta -i <archivo_instancia>\n";
        print_usage();
        return 1;
    }
    string archivo_instancia = argv[3];

    bool probabilista = (metodo == "Greedy-probabilista");
    float alpha = 0.0f;
    if (probabilista) {
        if (argc < 5) {
            cerr << "Error: falta <alpha>\n";
            print_usage();
            return 1;
        }
        alpha = stof(argv[4]);
        if (alpha <= 0.0f || alpha > 1.0f) {
            cerr << "Error: alpha debe estar en (0,1].\n";
            return 1;
        }
    } else if (metodo != "Greedy") {
        cerr << "Método desconocido: " << metodo << "\n";
        print_usage();
        return 1;
    }

    load_graph(archivo_instancia);
    if (n == 0) {
        cerr << "No se cargó el grafo.\n";
        return 1;
    }

    vector<int> solucion = probabilista
        ? misp_greedy_aleatorizado(alpha)
        : misp_heuristica_greedy();

    cout << "Solución (" << solucion.size() << " nodos):\n";
    for (int u : solucion) cout << u << " ";
    cout << "\n";
    return 0;
}