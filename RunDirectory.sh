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

# 2. Create results directory
RESULTS_DIR="$RUNDIR/results/$RUNNAME"
mkdir -p "$RESULTS_DIR"

# 3. Copy source files for record-keeping
cp src/Geometries/GDEntireDomainWithParticularities.* "$RESULTS_DIR/"

echo "=== Source files copied to $RESULTS_DIR ==="

# 4. Run ECOGEN with nohup in background
nohup mpirun -np $NP --use-hwthread-cpus "$ECOGEN_EXEC" "$RUNDIR" > "$RESULTS_DIR/output.log" 2>&1 &

SIM_PID=$!
echo "Simulation started in background (PID=$SIM_PID). Output is in $RESULTS_DIR/output.log"

# 5. Follow the output live (Ctrl-C to stop watching, simulation keeps running)
tail -f "$RESULTS_DIR/output.log"
