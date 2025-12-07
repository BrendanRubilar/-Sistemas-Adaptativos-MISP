#!/usr/bin/env bash
# Versión final y robusta del script runner.

# Termina el script si un comando falla (excepto donde se maneje explícitamente).
set -euo pipefail

# --- 1. Parámetros de Entrada ---
if [[ $# -lt 2 ]]; then
  echo "Uso: runner.sh <instancia> <semilla> [parámetros...]" >&2
  exit 1
fi
INSTANCE_FILE="$1"
SEED="$2"
shift 2
CONFIG_PARAMS="$@"

# --- 2. Configuración ---
TIME_LIMIT="10" # Tiempo para cada ejecución de irace
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXECUTABLE="$PROJECT_ROOT/MetaheuristicaPoblacional"

# --- 3. Ejecución ---
OUTPUT_FILE=$(mktemp)
trap 'rm -f "$OUTPUT_FILE"' EXIT

# Ejecuta el programa con timeout. El '|| true' es CRUCIAL para que
# 'set -e' no mate el script si el programa es terminado por timeout.
/usr/bin/timeout -k 2s "${TIME_LIMIT}s" \
    "$EXECUTABLE" BRKGA -i "$INSTANCE_FILE" -t "$TIME_LIMIT" -s "$SEED" $CONFIG_PARAMS \
    > "$OUTPUT_FILE" || true

# --- 4. Extracción del Resultado ---
# Se busca la calidad. El '|| true' al final es CRUCIAL para que 'set -e'
# no mate el script si 'grep' no encuentra ninguna solución.
FINAL_QUALITY=$(grep 'Calidad:' "$OUTPUT_FILE" | awk '{print $2}' | sed 's/,//' | tail -n 1 || true)

# Si no se encontró calidad (el programa falló o no encontró solución), se asigna 0.
if [[ -z "$FINAL_QUALITY" ]]; then
    FINAL_QUALITY=0
fi

# --- 5. Devolver el Costo a irace ---
echo "$(( -FINAL_QUALITY ))"