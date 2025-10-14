// simulation.js
import { plotEnergy, plotPositions3D, plotTemperature, plotPressure, plotMetaBias, plotProjection2D, plotCVHistogram, plotParticles3D } from "./plotlyUtils.js";
import { setButtonsDisabled } from "./ui.js";

let lastSimulationData = null;

/**
 * Establece un indicador de carga (spinner) en un elemento HTML dado.
 * @param {string} id - El ID del elemento donde mostrar el spinner.
 * @param {boolean} isLoading - True para mostrar el spinner, false para ocultarlo.
 *
 * Sets a loading indicator (spinner) on a given HTML element.
 * @param {string} id - The ID of the element where the spinner will be shown.
 * @param {boolean} isLoading - True to show the spinner, false to hide it.
 */
export const setLoading = (id, isLoading) => {
  const el = document.getElementById(id);
  if (!el) return;

  el.innerHTML = isLoading
    ? `<div class="flex justify-center items-center w-full h-full">
         <div class="w-10 h-10 border-4 border-t-transparent border-blue-500 rounded-full animate-spin"></div>
       </div>`
    : "";
};


/**
 * Obtiene el valor de un input HTML y lo convierte usando un parser.
 * Por defecto, convierte a número flotante.
 * @param {string} id - El ID del input HTML.
 * @param {Function} parser - Función para convertir el valor (por defecto parseFloat).
 * @returns {*} - El valor parseado del input.
 *
 * Retrieves the value of an HTML input and converts it using a parser.
 * By default, it converts the value to a floating-point number.
 * @param {string} id - The ID of the HTML input.
 * @param {Function} parser - Function to convert the value (default is parseFloat).
 * @returns {*} - The parsed value of the input.
 */
export const getInputValue = (id, parser = parseFloat) => parser(document.getElementById(id)?.value);

/**
 * Descarga los datos de la última simulación en un archivo JSON.
 * Muestra un alerta si no hay datos disponibles.
 *
 * @example
 * downloadSimulation(); // Descarga 'simulacion_2025-10-11T15-30-45-123Z.json'
 *
 * Downloads the last simulation data as a JSON file.
 * Shows an alert if no data is available.
 *
 * @example
 * downloadSimulation(); // Downloads 'simulacion_2025-10-11T15-30-45-123Z.json'
 */
export const downloadSimulation = () => {
  if (!lastSimulationData) return alert("⚠️ There are no simulations available to download. Please run one first.");

  const url = URL.createObjectURL(new Blob([JSON.stringify(lastSimulationData, null, 2)], { type: "application/json" }));

  Object.assign(document.createElement("a"), {
    href: url,
    download: `simulacion_${new Date().toISOString().replace(/[:.]/g, "-")}.json`,
  }).click();

  URL.revokeObjectURL(url);
};


/**
 * Genera el payload para enviar al backend según el método de simulación seleccionado.
 * Dependiendo del método, recoge los valores de los inputs correspondientes y los convierte
 * a los tipos adecuados (números enteros o flotantes).
 *
 * @param {string} method - Método de simulación seleccionado ('NVT_Langevin_BAOAB' o 'muVT_Monte_Carlo').
 * @returns {object} Payload listo para enviar al servidor con los parámetros de simulación.
 *
 * Generates the payload to send to the backend based on the selected simulation method.
 * Depending on the method, it collects input values and converts them to the appropriate types
 * (integers or floats).
 *
 * @param {string} method - Selected simulation method ('NVT_Langevin_BAOAB' or 'muVT_Monte_Carlo').
 * @returns {object} Payload ready to be sent to the server with the simulation parameters.
 */
export const createPayload = (method) =>
  method === "NVT_Langevin_BAOAB"
    ? {
        n_particles: getInputValue("n_particles", parseInt),
        n_steps: getInputValue("n_steps", parseInt),
        dt: getInputValue("dt"),
        box_size: getInputValue("box_size"),
        cutoff: getInputValue("cutoff"),
        skin: getInputValue("skin"),
        seed: getInputValue("seed", parseInt),
        kT: getInputValue("kT"),
        gamma: getInputValue("gamma"),
      }
    : method === "muVT_Monte_Carlo"
      ? (() => {
          const volume = getInputValue("V");
          return {
            n_particles: getInputValue("N", parseInt),
            volume,
            T: getInputValue("T"),
            mu: getInputValue("mu"),
            n_steps: getInputValue("nsteps", parseInt),
            xi_target: getInputValue("xi_target"),
            k_bias: getInputValue("k_bias"),
            sigma_meta: getInputValue("sigma_meta"),
            height_meta: getInputValue("height_meta"),
            box_size: Number.isFinite(volume) ? Math.cbrt(volume) : null,
          };
        })()
      : {};


/**
 * Genera los gráficos correspondientes a la simulación según el método seleccionado.
 * Dependiendo del método, actualiza los elementos HTML asociados y llama a las funciones de plot
 * correspondientes para mostrar los datos en gráficos 2D y 3D.
 *
 * @param {string} method - Método de simulación ('NVT_Langevin_BAOAB' o 'muVT_Monte_Carlo').
 * @param {object} data - Datos devueltos por la simulación para graficar.
 * @param {number} n_particles - Número de partículas en la simulación (usado en algunos plots).
 * @param {number} box_size - Tamaño de la caja de simulación (usado en algunos plots 3D).
 * @param {boolean} isDark - Indica si el tema es oscuro para ajustar los colores de los gráficos.
 *
 * Example:
 * plotSimulationData("NVT_Langevin_BAOAB", simulationData, 1000, 10.0, true);
 *
 * Plots the simulation data according to the selected method.
 * Depending on the method, it updates the associated HTML elements and calls the plotting functions
 * to render 2D and 3D charts.
 *
 * @param {string} method - Simulation method ('NVT_Langevin_BAOAB' or 'muVT_Monte_Carlo').
 * @param {object} data - Data returned from the simulation to be plotted.
 * @param {number} n_particles - Number of particles in the simulation (used in some plots).
 * @param {number} box_size - Simulation box size (used in some 3D plots).
 * @param {boolean} isDark - Indicates if dark theme is active to adjust chart colors.
 *
 * Example:
 * plotSimulationData("muVT_Monte_Carlo", simulationData, 500, 8.0, false);
 */
export const plotSimulationData = (method, data, n_particles, box_size, isDark) => {
  const config = {
    NVT_Langevin_BAOAB: {
      elements: ["output", "positions3d", "temperature", "pressure"],
      plots: [() => plotEnergy(data, isDark), () => plotPositions3D(data, n_particles, isDark), () => plotTemperature(data, n_particles, isDark), () => plotPressure(data, n_particles, box_size, isDark)],
    },
    muVT_Monte_Carlo: {
      elements: ["plotMetaBias", "plotParticles3D", "plotProjection2D", "plotCVHistogram"],
      plots: [() => plotMetaBias(data, isDark), () => plotParticles3D(data, isDark), () => plotProjection2D(data, isDark), () => plotCVHistogram(data, isDark)],
    },
  };

  const methodConfig = config[method];
  if (!methodConfig) return;

  methodConfig.elements.forEach((id) => setLoading(id, false));
  methodConfig.plots.forEach((fn) => fn());
};


/**
 * Ejecuta la simulación según el método seleccionado en la interfaz.
 *
 * 1. Obtiene el método de simulación del select HTML.
 * 2. Muestra los spinners de carga en todos los elementos de plot.
 * 3. Genera el payload con los parámetros de la simulación.
 * 4. Envía la solicitud POST al backend.
 * 5. Guarda los datos devueltos en `lastSimulationData`.
 * 6. Muestra errores si existen.
 * 7. Calcula número de partículas y tamaño de caja.
 * 8. Llama a `plotSimulationData` para graficar los resultados.
 * 9. Habilita los botones al finalizar.
 *
 * @async
 * @returns {Promise<void>} - No devuelve ningún valor.
 *
 * Example:
 * runSimulation(); // Ejecuta la simulación según los inputs actuales.
 *
 * Runs the simulation according to the selected method in the interface.
 *
 * 1. Gets the simulation method from the HTML select.
 * 2. Shows loading spinners on all plot elements.
 * 3. Creates the payload with simulation parameters.
 * 4. Sends a POST request to the backend.
 * 5. Stores the returned data in `lastSimulationData`.
 * 6. Displays errors if any.
 * 7. Calculates number of particles and box size.
 * 8. Calls `plotSimulationData` to render the results.
 * 9. Re-enables buttons when finished.
 *
 * @async
 * @returns {Promise<void>} - Does not return any value.
 *
 * Example:
 * runSimulation(); // Runs the simulation based on the current inputs.
 */
export const runSimulation = async () => {
  const method = document.getElementById("simulation_method")?.value;
  const messageDiv = document.getElementById("message");
  if (!messageDiv) return;
  messageDiv.innerHTML = "";

  ["output", "positions3d", "temperature", "pressure", "plotMetaBias", "plotParticles3D", "plotProjection2D", "plotCVHistogram"].forEach((id) => setLoading(id, true));

  const payload = createPayload(method);
  setButtonsDisabled(true);

  try {
    const res = await fetch(`/simulate/${method}`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });

    const data = await res.json();
    lastSimulationData = data;

    if (data.error) {
      messageDiv.innerHTML = `<p style="color:red;">${data.error}</p>`;
      return;
    }

    const n_particles = data.positions?.[0]?.[0]?.length ?? payload.n_particles;
    const box_size = payload.box_size ?? Math.cbrt(payload.volume ?? 0) ?? data.box_size;
    const isDark = document.documentElement.classList.contains("dark");

    plotSimulationData(method, data, n_particles, box_size, isDark);
  } catch (err) {
    console.error(err);
    messageDiv.innerHTML = `<p style="color:red;">Error: ${err?.message}</p>`;
  } finally {
    setButtonsDisabled(false);
  }
};


/**
 * Gestiona el cambio de método de simulación en la interfaz.
 * Muestra u oculta los inputs correspondientes a cada tipo de simulación (NVT o muVT)
 * según la opción seleccionada en el select.
 *
 * @example
 * // Selecciona 'muVT_Monte_Carlo', se muestran los inputs de muVT y se ocultan los de NVT
 * simulationSelect.dispatchEvent(new Event('change'));
 *
 * Manages the simulation method change in the interface.
 * Shows or hides the inputs corresponding to each simulation type (NVT or muVT)
 * depending on the selected option in the select element.
 *
 * @example
 * // Selects 'muVT_Monte_Carlo', shows muVT inputs, hides NVT inputs
 * simulationSelect.dispatchEvent(new Event('change'));
 */
const simulationSelect = document.getElementById("simulation_method");
const muvDivs = document.querySelectorAll(".simulation-muvt");
const nvtDivs = document.querySelectorAll(".simulation-nvt");
const nvtInputs = document.getElementById("inputs_nvt");
const muvInputs = document.getElementById("inputs_muvt");

simulationSelect.addEventListener("change", () => {
  const isMUVT = simulationSelect.value === "muVT_Monte_Carlo";

  [...muvDivs, muvInputs].forEach((el) => el.classList.toggle("hidden", !isMUVT));
  [...nvtDivs, nvtInputs].forEach((el) => el.classList.toggle("hidden", isMUVT));
});
