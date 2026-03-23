#!/bin/bash

# --- CONFIGURAÇÃO ---
# Lista dos nós (Exemplo: ap20n0-ib0, ap20n0-ib1, etc.)
NOS=("ap20n0-ib0" "ap20n1-ib0" "ap20n2-ib0" "ap20n3-ib0")

# Caminhos
PYTHON_EXE="/home/gilberto/miniconda3/bin/python"
SCRIPT="/scratch/gilberto/scripts/indices_icon.py"

TOTAL_NOS=${#NOS[@]}

echo "Iniciando Cluster ICON: $TOTAL_NOS nós detectados."

for i in "${!NOS[@]}"; do
    NODE=${NOS[$i]}
    RANK=$i
    
    echo " -> Lançando Job no nó: $NODE (Rank $RANK)..."
    
    # Comando SSH:
    # 1. Entra no nó
    # 2. Roda nohup (background)
    # 3. Passa --rank e --size para o python saber qual parte processar
    ssh $NODE "nohup $PYTHON_EXE $SCRIPT --rank $RANK --size $TOTAL_NOS > /scratch/gilberto/logs/log_icon_no${RANK}.txt 2>&1 &"

done

echo "---------------------------------------------------"
echo "Sucesso! Logs disponíveis em /scratch/gilberto/logs/log_icon_noX.txt"
