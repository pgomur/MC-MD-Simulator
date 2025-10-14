import { initUI } from "./ui.js";
import { runSimulation, downloadSimulation } from "./simulation.js";

/**
 * Inicializa la interfaz de usuario y asigna los eventos a los botones
 * cuando el DOM ha sido completamente cargado.
 *
 * Funcionalidades:
 * - Llama a `initUI()` para configurar toda la UI.
 * - Asocia los botones "runSimulation" y "downloadSimulation" a sus funciones correspondientes.
 *
 * @example
 * // Se ejecuta automáticamente al cargar la página
 * document.addEventListener("DOMContentLoaded", ...);
 *
 * Initializes the user interface and attaches button events
 * when the DOM is fully loaded.
 *
 * Features:
 * - Calls `initUI()` to set up the entire UI.
 * - Attaches the "runSimulation" and "downloadSimulation" buttons to their respective functions.
 *
 * @example
 * // Automatically runs when the page is loaded
 * document.addEventListener("DOMContentLoaded", ...);
 */
document.addEventListener("DOMContentLoaded", () => {
  initUI();

  const buttonActions = {
    runSimulation,
    downloadSimulation,
  };

  for (const [id, action] of Object.entries(buttonActions)) {
    const btn = document.getElementById(id);
    if (btn) btn.addEventListener("click", action);
  }
});
