#!/bin/bash


CORES=1

DIRECTORIO_BASE=/home/microbios/silvia_backup/silvia/micro/2021_06_28__null_models/results_null_model/2024_simcomms_WITH_GROUPS

# Verificar que el directorio existe
if [ ! -d "$DIRECTORIO_BASE" ]; then
    echo "El directorio $DIRECTORIO_BASE no existe."
    exit 1
fi

# Iterar sobre todas las carpetas dentro del directorio base
for carpeta in "$DIRECTORIO_BASE"/SIM*; do
    if [ -d $carpeta ]; then
        echo "Procesando carpeta: $carpeta"

	# Get the absolute path to the reorganizar script
	SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/02_process_results"
	
	echo '#!/bin/bash' > launch_$(basename $carpeta).sh
	echo 'python '$SCRIPT_DIR'/reorganizar_comunidades.py '$carpeta >> launch_$(basename $carpeta).sh
	chmod +x launch_$(basename $carpeta).sh

	qsubmit.pl -q x86_64 -n $CORES --mem 20 -s launch_$(basename $carpeta).sh
    fi
done

echo "Procesamiento lanzado para todas las carpetas en $DIRECTORIO_BASE."

