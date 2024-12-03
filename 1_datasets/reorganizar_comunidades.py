import os
import pandas as pd
import re

# Configuración de argumentos de línea de comandos
parser = argparse.ArgumentParser(description="Reorganiza los archivos CSV de comunidades por puntos temporales.")
parser.add_argument(
    "source_dir",
    help="Directorio que contiene los archivos CSV originales."
)
args = parser.parse_args()

# Variables de directorio
source_dir = args.source_dir
dest_dir = "new_" + os.path.basename(source_dir)

# Crear la carpeta de destino si no existe
if not os.path.exists(dest_dir):
    os.makedirs(dest_dir)

# Función para extraer el número <n> del nombre del archivo
def extract_number(filename):
    match = re.search(r"simul_1_t_(\d+)\.csv", filename)
    return int(match.group(1)) if match else float('inf')

# Obtener todos los archivos del directorio fuente y ordenarlos numéricamente
file_list = sorted(
    [f for f in os.listdir(source_dir) if f.startswith("simul_1_t_") and f.endswith(".csv")],
    key=extract_number
)

# Leer todos los archivos y almacenarlos en una lista de DataFrames
dataframes = []
for file in file_list:
    file_path = os.path.join(source_dir, file)
    df = pd.read_csv(file_path, usecols=lambda column: column != df.columns[0]) # drop first column!
    dataframes.append(df)

# Número de comunidades (filas en cada archivo, excluyendo la cabecera)
num_communities = len(dataframes[0])  # Se asume que todos los archivos tienen las mismas filas

# Reorganizar datos: crear un archivo por comunidad
for community_idx in range(num_communities):
    community_data = []
    for df in dataframes:
        community_data.append(df.iloc[community_idx])  # Tomar la fila correspondiente de cada punto temporal

    # Combinar datos y modificar la primera columna
    combined_df = pd.DataFrame(community_data)
    combined_df.reset_index(drop=True, inplace=True)  # Reiniciar los índices
    combined_df.insert(0, "time", combined_df.index)  # Insertar la columna "time" con índices
    combined_df.columns.values[0] = "time"  # Asegurarse de que el nombre de la columna sea "time"

    # Guardar como nuevo archivo CSV
    output_file = os.path.join(dest_dir, f"community_{community_idx + 1}.csv")
    combined_df.to_csv(output_file, index=False)

print(f"Reorganización completada. Archivos creados en: {dest_dir}")

