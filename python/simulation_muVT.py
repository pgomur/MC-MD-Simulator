# simulation_muVT.py
import os
import subprocess
import struct
import numpy as np
import json
from flask import jsonify

def validate_params_muVT(params):
    """
    Valida los parámetros de la simulación μVT y aplica valores por defecto.
    Validate μVT simulation parameters and apply default values.
    """
    try:
        return {
            "n_particles": max(1, int(params.get("n_particles", 20))),      # Número de partículas / Number of particles
            "V": max(1.0, float(params.get("V", 300.0))),                  # Volumen del sistema / System volume
            "T": max(0.01, float(params.get("T", 1.0))),                   # Temperatura / Temperature
            "mu": float(params.get("mu", 0.5)),                             # Potencial químico / Chemical potential
            "n_steps": max(1, int(params.get("n_steps", 5000))),           # Número de pasos MC / Number of MC steps
            "xi_target": float(params.get("xi_target", 0.5)),              # Umbrella sampling objetivo / Umbrella sampling target
            "k_bias": float(params.get("k_bias", 10.0)),                   # Constante de bias / Bias constant
            "sigma_meta": float(params.get("sigma_meta", 0.1)),            # Parámetro metadynamics / Metadynamics sigma
            "height_meta": float(params.get("height_meta", 1.0))           # Altura del bias metadynamics / Metadynamics height
        }
    except (ValueError, TypeError) as e:
        raise ValueError(f"Parámetros inválidos: {str(e)} / Invalid parameters: {str(e)}")

def read_binary_muVT(file_path):
    """
    Lee el archivo binario de la simulación μVT Monte Carlo y devuelve un diccionario de datos.
    Read μVT Monte Carlo binary file and return a dictionary of data.
    """
    try:
        with open(file_path, "rb") as f:
            # 1. Número de partículas / Number of particles
            n_particles = struct.unpack('i', f.read(4))[0]

            # 2-4. energy, density, max_disp
            energy = struct.unpack('d', f.read(8))[0]
            density = struct.unpack('d', f.read(8))[0]
            max_disp = struct.unpack('d', f.read(8))[0]

            # 5. Umbrella sampling
            umbrella_on = struct.unpack('i', f.read(4))[0] != 0
            xi_target_bin = k_bias_bin = None
            if umbrella_on:
                xi_target_bin = struct.unpack('d', f.read(8))[0]
                k_bias_bin = struct.unpack('d', f.read(8))[0]

            # 6. Metadynamics
            meta_on = struct.unpack('i', f.read(4))[0] != 0
            sigma_meta_bin = height_meta_bin = None
            meta_bias = None
            if meta_on:
                sigma_meta_bin = struct.unpack('d', f.read(8))[0]
                height_meta_bin = struct.unpack('d', f.read(8))[0]
                meta_bias = np.frombuffer(f.read(n_particles*8), dtype='d')

            # Posiciones de partículas / Particle positions
            positions = np.frombuffer(f.read(3*n_particles*8), dtype='d').reshape((3, n_particles))

    except (struct.error, IOError, OSError) as e:
        raise RuntimeError(f"Error leyendo binario μVT: {str(e)} / Error reading μVT binary: {str(e)}")

    return {
        "n_particles": n_particles,
        "energy": energy,
        "density": density,
        "max_disp": max_disp,
        "umbrella_on": umbrella_on,
        "xi_target": xi_target_bin,
        "k_bias": k_bias_bin,
        "meta_on": meta_on,
        "sigma_meta": sigma_meta_bin,
        "height_meta": height_meta_bin,
        "meta_bias": meta_bias.tolist() if meta_bias is not None else None,
        "positions": positions.tolist()
    }

def run_simulation_muVT(params, BASE_DIR):
    """
    Ejecuta la simulación μVT Monte Carlo y devuelve los resultados en JSON.
    Run μVT Monte Carlo simulation and return results as JSON.
    """
    # Validar parámetros / Validate parameters
    try:
        p = validate_params_muVT(params)
    except ValueError as e:
        return jsonify({"error": str(e)}), 400

    # Rutas de ejecutable y archivos / Paths for executable and output files
    FORTRAN_EXE = os.path.join(BASE_DIR, "fortran", "MuVT_MonteCarlo", "bin", "mc_muvt")
    BINARY_PATH = os.path.join(BASE_DIR, "data", "MuVT_MonteCarlo", "simulation_muVT.bin")
    PREDICTION_JSON = os.path.join(BASE_DIR, "data", "MuVT_MonteCarlo", "predictions_muVT.json")

    if not os.path.exists(FORTRAN_EXE):
        return jsonify({"error": f"Ejecutable no encontrado: {FORTRAN_EXE} / Executable not found: {FORTRAN_EXE}"}), 500

    # Crear carpeta si no existe / Create output folder if missing
    os.makedirs(os.path.dirname(BINARY_PATH), exist_ok=True)

    # Construir comando de ejecución / Build execution command
    cmd = [
        FORTRAN_EXE,
        str(p["n_particles"]),
        str(p["V"]),
        str(p["T"]),
        str(p["mu"]),
        str(p["n_steps"]),
        BINARY_PATH,
        str(p["xi_target"]),
        str(p["k_bias"]),
        str(p["sigma_meta"]),
        str(p["height_meta"])
    ]
    print("Ejecutando comando μVT Monte Carlo / Running μVT Monte Carlo command:", " ".join(cmd), flush=True)

    # Ejecutar simulación / Run simulation
    try:
        result = subprocess.run(
            cmd,
            cwd=os.path.dirname(BINARY_PATH),
            check=True,
            capture_output=True,
            text=True
        )
        if result.stderr:
            print("Fortran stderr:", result.stderr, flush=True)
    except subprocess.CalledProcessError as e:
        return jsonify({"error": f"Error en simulación μVT: {e.stderr} / μVT simulation error: {e.stderr}"}), 500

    if not os.path.exists(BINARY_PATH):
        return jsonify({"error": f"Archivo binario no creado: {BINARY_PATH} / Binary file not created: {BINARY_PATH}"}), 500

    # Leer binario / Read binary file
    try:
        data_json = read_binary_muVT(BINARY_PATH)
    except RuntimeError as e:
        return jsonify({"error": str(e)}), 500

    # Guardar JSON / Save JSON
    with open(PREDICTION_JSON, "w") as f:
        json.dump(data_json, f, indent=2)

    return jsonify(data_json)
