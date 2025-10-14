from flask import Flask, request, jsonify
from pathlib import Path
from simulation_NVT import run_simulation
from simulation_muVT import run_simulation_muVT

# Definir la ruta base del proyecto usando pathlib
# Define the base directory of the project using pathlib
BASE_DIR = Path(__file__).resolve().parent.parent

# Definir la carpeta donde están los archivos web estáticos (index.html, etc.)
# Define the folder containing static web files (index.html, etc.)
WEB_DIR = BASE_DIR / "web"

# Crear la aplicación Flask, indicando la carpeta de archivos estáticos
# Create the Flask app, specifying the static files folder
app = Flask(__name__, static_folder=str(WEB_DIR), static_url_path="")

# Ruta principal "/" que devuelve el index.html
# Main route "/" that returns the index.html
@app.route("/")
def index():
    return app.send_static_file("index.html")

# Ruta para ejecutar simulaciones mediante POST
# Route to execute simulations via POST
@app.route("/simulate/<simulation_type>", methods=["POST"])
def simulate(simulation_type: str):
    """
    Maneja solicitudes de simulación usando pattern-matching y desestructuración de diccionarios.
    Handle simulation requests using pattern matching and dictionary destructuring.
    """
    # Obtener datos JSON del request
    # Get JSON data from the request
    params = request.get_json()

    # Validar que el JSON es un diccionario
    # Validate that the JSON is a dictionary
    if not isinstance(params, dict):
        return jsonify({"error": "Missing or invalid JSON payload"}), 400

    # Selección de la simulación usando match/case
    # Select the simulation using match/case
    match simulation_type, params:
        # Simulación NVT Langevin con algoritmo BAOAB
        # NVT Langevin simulation using BAOAB algorithm
        case "NVT_Langevin_BAOAB", p if isinstance(p, dict):
            return run_simulation(p, BASE_DIR)

        # Simulación muVT Monte Carlo
        # muVT Monte Carlo simulation
        case "muVT_Monte_Carlo", p if isinstance(p, dict):
            return run_simulation_muVT(p, BASE_DIR)

        # Si el tipo de simulación no está soportado
        # If the simulation type is not supported
        case _:
            return jsonify({"error": "Simulation type not supported"}), 400

# Ejecutar la aplicación si este archivo se ejecuta directamente
# Run the app if this file is executed directly
if __name__ == "__main__":
    # debug=True solo para desarrollo / debug=True only for development
    app.run(host="0.0.0.0", port=5000, debug=True)
