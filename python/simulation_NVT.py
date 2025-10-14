# simulation.py
import os
import subprocess
import struct
import numpy as np
import json
from flask import jsonify

def validate_params(params):
    """
    Valida los parámetros de simulación y aplica valores por defecto.
    Validate simulation parameters and apply default values.
    """
    try:
        return {
            "n_particles": max(1, int(params.get("n_particles", 10000))),  # Número de partículas / Number of particles
            "n_steps": max(1, int(params.get("n_steps", 2000))),           # Número de pasos de simulación / Number of simulation steps
            "dt": max(0.0001, float(params.get("dt", 0.005))),            # Paso de tiempo / Time step
            "box_size": max(1.0, float(params.get("box_size", 50.0))),    # Tamaño de la caja / Box size
            "cutoff": max(0.1, float(params.get("cutoff", 2.5))),         # Distancia de corte / Cutoff distance
            "skin": max(0.01, float(params.get("skin", 0.3))),            # Parámetro skin / Skin parameter
            "seed": max(1, int(params.get("seed", 12345))),               # Semilla aleatoria / Random seed
            "kT": float(params.get("kT", 1.0)),                           # Temperatura / Temperature
            "gamma": float(params.get("gamma", 0.1))                      # Coeficiente de fricción / Friction coefficient
        }
    except (ValueError, TypeError) as e:
        raise ValueError(f"Parámetros inválidos: {str(e)}")  # Error en validación / Validation error

def read_binary(file_path, n_steps, n_particles):
    """
    Lee el archivo binario de la simulación y devuelve una lista de frames.
    Read the simulation binary file and return a list of frames.
    """
    frames = []

    try:
        with open(file_path, "rb") as f:
            for step in range(n_steps):
                # Leer marcador de inicio de registro (4 bytes) / Read start marker (4 bytes)
                marker_start = f.read(4)
                if len(marker_start) != 4:
                    break
                record_size = struct.unpack('i', marker_start)[0]

                # Leer datos del registro / Read record data
                data = f.read(record_size)
                if len(data) != record_size:
                    break

                # Leer posiciones, velocidades y fuerzas / Read positions, velocities, and forces
                x = np.frombuffer(data[0:24*n_particles], dtype='d').reshape((3, n_particles))
                v = np.frombuffer(data[24*n_particles:48*n_particles], dtype='d').reshape((3, n_particles))
                f_array = np.frombuffer(data[48*n_particles:72*n_particles], dtype='d').reshape((3, n_particles))

                # Leer energías cinética, potencial y total / Read kinetic, potential, and total energies
                ke = struct.unpack('d', data[72*n_particles:72*n_particles+8])[0]
                pe = struct.unpack('d', data[72*n_particles+8:72*n_particles+16])[0]
                total = struct.unpack('d', data[72*n_particles+16:72*n_particles+24])[0]

                # Leer marcador de fin de registro / Read end marker
                marker_end = f.read(4)
                if len(marker_end) != 4:
                    break
                if struct.unpack('i', marker_end)[0] != record_size:
                    raise ValueError(f"Marcadores inconsistentes en paso {step+1} / Inconsistent markers at step {step+1}")

                # Agregar frame a la lista / Append frame to list
                frames.append({
                    'step': step + 1,
                    'total_energy': float(total),
                    'kinetic_energy': float(ke),
                    'potential_energy': float(pe),
                    'positions': x.tolist(),
                    'forces': f_array.tolist()
                })
    except (struct.error, IOError, OSError) as e:
        raise RuntimeError(f"Error leyendo binario: {str(e)} / Error reading binary: {str(e)}")

    if len(frames) != n_steps:
        raise RuntimeError(f"Número de frames incorrecto. Esperado: {n_steps}, Leídos: {len(frames)} / Incorrect number of frames. Expected: {n_steps}, Read: {len(frames)}")

    return frames

def run_simulation(params, BASE_DIR):
    """
    Ejecuta la simulación y devuelve resultados como JSON.
    Run the simulation and return results as JSON.
    """
    # Validar parámetros / Validate parameters
    try:
        p = validate_params(params)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400

    # Rutas de ejecutable y archivos / Paths for executable and files
    FORTRAN_EXE = os.path.join(BASE_DIR, "fortran", "NVT_Langevin_BAOAB", "bin", "md_baoab")
    BINARY_PATH = os.path.join(BASE_DIR, "data", "NVT_Langevin_BAOAB", "simulation_NVT_BAOAB.bin")
    PREDICTION_JSON = os.path.join(BASE_DIR, "data", "NVT_Langevin_BAOAB", "predictions_NVT_BAOAB.json")

    if not os.path.exists(FORTRAN_EXE):
        return jsonify({"error": f"Ejecutable no encontrado: {FORTRAN_EXE} / Executable not found: {FORTRAN_EXE}"}), 500

    # Preparar directorios y borrar archivos antiguos / Prepare directories and remove old files
    os.makedirs(os.path.dirname(BINARY_PATH), exist_ok=True)
    for file_path in [BINARY_PATH, PREDICTION_JSON]:
        if os.path.exists(file_path):
            os.remove(file_path)

    # Construir comando / Build command
    cmd = [
        FORTRAN_EXE,
        str(p["n_particles"]),
        str(p["n_steps"]),
        str(p["dt"]),
        str(p["box_size"]),
        str(p["cutoff"]),
        str(p["skin"]),
        str(p["seed"]),
        str(p["kT"]),
        str(p["gamma"])
    ]
    print("Ejecutando comando simulate / Running simulation command:", " ".join(cmd), flush=True)

    # Ejecutar simulación / Run simulation
    try:
        result = subprocess.run(cmd, cwd=os.path.dirname(BINARY_PATH), check=True,
                                capture_output=True, text=True)
        if result.stderr:
            print(f"Fortran stderr: {result.stderr}", flush=True)
    except subprocess.CalledProcessError as e:
        return jsonify({"error": f"Error en simulación: {e.stderr} / Simulation error: {e.stderr}"}), 500

    # Verificar que el binario fue creado / Check binary file created
    if not os.path.exists(BINARY_PATH):
        return jsonify({"error": f"Archivo binario no creado: {BINARY_PATH} / Binary file not created: {BINARY_PATH}"}), 500

    # Leer binario / Read binary file
    try:
        frames = read_binary(BINARY_PATH, p["n_steps"], p["n_particles"])
    except (RuntimeError, ValueError) as e:
        return jsonify({"error": str(e)}), 500

    # Preparar datos de salida / Prepare output data
    output_data = {
        "steps": [f['step'] for f in frames],
        "total_energy": [f['total_energy'] for f in frames],
        "kinetic_energy": [f['kinetic_energy'] for f in frames],
        "potential_energy": [f['potential_energy'] for f in frames],
        "positions": [f['positions'] for f in frames],
        "forces": [f['forces'] for f in frames]
    }

    # Guardar JSON / Save JSON
    with open(PREDICTION_JSON, "w") as f:
        json.dump(output_data, f, indent=2)

    # Retornar resultado JSON para Flask / Return JSON response for Flask
    return jsonify({
        "steps": output_data["steps"],
        "labels": output_data["total_energy"],
        **output_data
    })
