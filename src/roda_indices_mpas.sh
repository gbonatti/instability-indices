#!/bin/bash

# --- CONFIGURAÇÃO ---
# Lista dos nós do Cluster
NOS=("ap20n4-ib0" "ap20n5-ib0" "ap20n6-ib0" "ap20n7-ib0")

# Caminhos
PYTHON_EXE="/home/gilberto/miniconda3/bin/python"
SCRIPT="/scratch/gilberto/scripts/indices_mpas.py"

TOTAL_NOS=${#NOS[@]}

echo "Iniciando Processamento MPAS Distribuído ($TOTAL_NOS nós)..."

for i in "${!NOS[@]}"; do
    NODE=${NOS[$i]}
    RANK=$i
    
    echo " -> Lançando Job no nó: $NODE (Rank $RANK)..."
    
    # Comando SSH:
    # O Python vai receber --rank e --size e saberá quais horas processar
    ssh $NODE "nohup $PYTHON_EXE $SCRIPT --rank $RANK --size $TOTAL_NOS > /scratch/gilberto/logs/log_mpas_no${RANK}.txt 2>&1 &"

done

echo "---------------------------------------------------"
echo "Jobs submetidos! Acompanhe: tail -f /scratch/gilberto/logs/log_mpas_no*.txt"
