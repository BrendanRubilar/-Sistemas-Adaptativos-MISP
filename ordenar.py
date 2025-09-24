import csv
import sys
from collections import defaultdict
from statistics import mean

def main():
    if len(sys.argv) < 2:
        print("Uso: python3 ordenar.py <input.csv> [output.csv]")
        sys.exit(1)

    in_path = sys.argv[1]
    out_path = sys.argv[2] if len(sys.argv) >= 3 else "resumen_densidades.csv"

    # Acumular por densidad
    by_density_means = defaultdict(list)    # densidad -> lista de t_mean
    by_density_stdevs = defaultdict(list)   # densidad -> lista de t_stdev

    with open(in_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            density = row["density"].strip()  # conservar representación textual
            try:
                t_mean = float(row["t_mean"])
                t_stdev = float(row["t_stdev"])
            except (KeyError, ValueError):
                continue
            by_density_means[density].append(t_mean)
            by_density_stdevs[density].append(t_stdev)

    # Escribir resumen
    # Ordenar densidades numéricamente
    densities_sorted = sorted(by_density_means.keys(), key=lambda d: float(d))

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["density", "t_mean_avg", "t_stdev_avg"])
        for d in densities_sorted:
            if by_density_means[d]:
                avg_mean = mean(by_density_means[d])
                avg_stdev = mean(by_density_stdevs[d])
                # Formato en notación científica similar al input
                writer.writerow([d, f"{avg_mean:.6e}", f"{avg_stdev:.6e}"])

    print(f"Resumen escrito en: {out_path}")

if __name__ == "__main__":
    main()