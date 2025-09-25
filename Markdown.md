# Heurísticas para MISP (Maximum Independent Set Problem)

## 1. Archivos principales
- `Heuristica.cpp`: Implementa heurística Greedy y Greedy-probabilista.
- `uhr.cpp`: Ejecuta experimentos sobre múltiples instancias (tamaños y densidades), midiendo tiempo y tamaño del conjunto independiente, se usa para el dataset completo.
- `dataset_grafos_no_dirigidos/`: Dataset a utilizar para las pruebas.

## 2. Requisitos
- Compilador: g++ >= C++17
- Sistema: Linux

## 3. Compilación heurísticas
```bash
g++ -std=c++17 -O2 -o Heuristica Heuristica.cpp
```

## 4. Uso del ejecutable
Formato general:
```
./Heuristica Greedy -i <ruta_instancia>
./Heuristica Greedy-probabilista -i <ruta_instancia> <alpha>
```
Parámetros:
- `<ruta_instancia>`: archivo .graph (primer número = n, luego pares u v).
- `<alpha>`: 0 < alpha ≤ 1. Controla el grado de aleatoriedad (porción de candidatos considerados).

Ejemplos:
```bash
./Heuristica Greedy -i erdos_n1000_p0c0.1_1.graph
./Heuristica Greedy-probabilista -i erdos_n1000_p0c0.1_1.graph 0.1
```

Salida:
```
Solución (k nodos):
<lista de nodos> 
```

## 5. Experimentos de dataset completo (Usar uhr.cpp)
Compilar:
```bash
g++ -std=c++17 -O2 -o uhr uhr.cpp
```
Uso:
```
./uhr <csv_salida> <runs> <instancia_inicial> <instancia_final> <paso>
```
Ejemplo (30 instancias, runs=32):
```bash
./uhr resultados.csv 32 1 30 1
```

Genera:
- `resultados.csv`: Resumen global por (size,density).
- `resultados_1000.csv`, `resultados_2000.csv`, `resultados_3000.csv`: Resumen por tamaño.

Columnas:
```
size,density,t_mean,t_stdev,iset_mean,iset_stdev
```
Tiempos en nanosegundos (std::chrono, high_resolution_clock).
