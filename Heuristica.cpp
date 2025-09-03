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
        cerr << "Error abriendo archivo\n";
        exit(1);
    }

    f >> n;
    graph.assign(n, Vertice{});

    int u, v;
    while (f >> u >> v) {
        graph[u].vecinos.push_back(v);
        graph[v].vecinos.push_back(u);
    }

    f.close();
}

int main() {
    load_graph("grafo_1.graph");

    for (int i = 0; i < 10; i++) {
        cout << "Nodo " << i << ": ";
        for (int v : graph[i].vecinos) cout << v << " ";
        cout << "\n";
    }

    return 0;
}

/*El plan es usar el algoritmo de Vertex Cover y usar su complemento para resolver esta shiet wola */