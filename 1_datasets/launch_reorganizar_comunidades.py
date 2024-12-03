#!/bin/bash

# Verifica si se proporcionó un argumento
if [ "$#" -ne 1 ]; then
    echo "Uso: $0 <directorio_base>"
    exit 1
fi

# Directorio base pasado como argumento
DIRECTORIO_BASE=$1

# Verificar que el directorio existe
if [ ! -d "$DIRECTORIO_BASE" ]; then
    echo "El directorio $DIRECTORIO_BASE no existe."
    exit 1
fi

# Iterar sobre todas las carpetas dentro del directorio base
for carpeta in "$DIRECTORIO_BASE"/*/; do
    if [ -d "$carpeta" ]; then
        echo "Procesando carpeta: $carpeta"
        # Llama al script de Python con la carpeta como argumento
        python reorganizar_comunidades.py "$carpeta"
    fi
done

echo "Procesamiento completado para todas las carpetas en $DIRECTORIO_BASE."

