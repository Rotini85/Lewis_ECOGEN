#!/bin/bash
# run_rt.sh — Run Rayleigh–Taylor ECOGEN case from rundirectory

# --- Config ---
RUNDIR="./rundirectory/RT/"
ECOGEN_EXEC="ECOGEN"
NP=14  # Number of MPI processes

# 1. Extract run name from main.xml
RUNNAME=$(grep -oP '(?<=<run>).*?(?=</run>)' "$RUNDIR/libtest/main.xml")

if [ -z "$RUNNAME" ]; then
    echo "Error: Could not extract <run> name from $RUNDIR/libtest/main.xml"
    exit 1
fi

echo "=== Starting ECOGEN simulation: $RUNNAME ==="

# 2. Create results directory (inside ECOGEN)
RESULTS_DIR="$RUNDIR/results/$RUNNAME"
mkdir -p "$RESULTS_DIR"

# 3. Create SimulationLogs directory OUTSIDE ECOGEN
PARENT_DIR="$(dirname "$(pwd)")"
SIMLOG_DIR="$PARENT_DIR/SimulationLogs/$RUNNAME"
mkdir -p "$SIMLOG_DIR"

# 4. Copy source files to SimulationLogs
cp src/Geometries/GDEntireDomainWithParticularities.* "$SIMLOG_DIR/"

echo "=== Source files copied to $SIMLOG_DIR ==="

# 5. Run ECOGEN with nohup in background
nohup mpirun -np $NP --use-hwthread-cpus "$ECOGEN_EXEC" "$RUNDIR" > "$RESULTS_DIR/output.log" 2>&1 &

SIM_PID=$!
echo "Simulation started in background (PID=$SIM_PID). Output is in $RESULTS_DIR/output.log"

# 6. Follow the output live (Ctrl-C to stop watching, simulation keeps running)
tail -f "$RESULTS_DIR/output.log"