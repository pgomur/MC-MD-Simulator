import { updatePlotlyTheme } from "./plotlyUtils.js";

/**
 * Habilita o deshabilita los botones de la interfaz y aplica clases de estilo
 * para reflejar el estado visualmente.
 *
 * @param {boolean} isDisabled - True para deshabilitar los botones, false para habilitarlos.
 *
 * @example
 * setButtonsDisabled(true); // Deshabilita los botones y aplica estilo de desactivado
 *
 * Enables or disables buttons in the interface and applies CSS classes
 * to visually reflect their state.
 *
 * @param {boolean} isDisabled - True to disable the buttons, false to enable them.
 *
 * @example
 * setButtonsDisabled(false); // Enables the buttons and restores normal style
 */
export function setButtonsDisabled(isDisabled) {
  ["runSimulation", "downloadSimulation"].forEach((id) => {
    const btn = document.getElementById(id);
    if (!btn) return;
    btn.disabled = isDisabled;
    btn.classList.toggle("opacity-50", isDisabled);
    btn.classList.toggle("cursor-not-allowed", isDisabled);
  });
}


/**
 * Carga de manera asíncrona el contenido HTML de un modal desde un archivo externo
 * y lo inserta en el contenedor del modal.
 *
 * @async
 * @returns {Promise<void>} - No devuelve ningún valor.
 *
 * @example
 * await loadModalHTML(); // Carga el modal en el contenedor
 *
 * Asynchronously loads the HTML content of a modal from an external file
 * and inserts it into the modal container.
 *
 * @async
 * @returns {Promise<void>} - Does not return any value.
 *
 * @example
 * await loadModalHTML(); // Loads the modal into the container
 */
export async function loadModalHTML() {
  const html = await (await fetch("./infoModal.html")).text();
  document.getElementById("modalContainer").innerHTML = html;
}


/**
 * Configura un overlay que bloquea la web en dispositivos móviles o pantallas pequeñas.
 *
 * - Detecta el tipo de dispositivo usando UAParser.
 * - Oculta el body y muestra un overlay si el dispositivo es móvil, tablet o la ventana es pequeña.
 * - El overlay informa al usuario que la web solo está disponible en escritorio.
 * - El overlay se elimina automáticamente si el tamaño de la ventana es suficiente.
 * - Se actualiza al cargar la página y al redimensionar la ventana.
 *
 * @example
 * setupDeviceOverlay(); // Activa el overlay para móviles y pantallas pequeñas
 *
 * Sets up an overlay that blocks the website on mobile devices or small screens.
 *
 * - Detects the device type using UAParser.
 * - Hides the body and shows an overlay if the device is mobile, tablet, or window is small.
 * - The overlay informs the user that the site is only available on desktop.
 * - The overlay is automatically removed if the window is large enough.
 * - Updates on page load and window resize events.
 *
 * @example
 * setupDeviceOverlay(); // Activates the overlay for mobile and small screens
 */
function setupDeviceOverlay() {
  const { device } = new UAParser().getResult();

  const checkDeviceAndWidth = () => {
    const body = document.body;
    const isMobileOrSmall = ["mobile", "tablet"].includes(device.type) || window.innerWidth <= 1023;

    const overlay = document.getElementById("mobileOverlay");

    if (isMobileOrSmall) {
      body.style.display = "none";

      if (!overlay) {
        const div = document.createElement("div");
        div.id = "mobileOverlay";
        div.className = "fixed top-0 left-0 w-screen h-screen bg-gray-900 text-white flex flex-col items-center justify-center p-6 z-[9999] text-center";
        div.innerHTML = `
          <h1 class="text-2xl font-bold mb-4">This website is not available on mobile or small screens.</h1>
          <p>Please use a desktop computer or increase the window size.</p>
        `;
        document.documentElement.appendChild(div);
      }
    } else {
      body.style.display = "";
      overlay?.remove();
    }
  };

  ["load", "resize"].forEach((evt) => window.addEventListener(evt, checkDeviceAndWidth));
}


/**
 * Configura el cambio de método de simulación en la interfaz.
 *
 * - Muestra u oculta los inputs correspondientes a cada método de simulación.
 * - Limpia los gráficos existentes para evitar confusión al cambiar de método.
 *
 * @example
 * setupSimulationMethodToggle(); // Activa el toggle de inputs y limpieza de gráficos
 *
 * Sets up the simulation method toggle in the interface.
 *
 * - Shows or hides the inputs corresponding to each simulation method.
 * - Clears existing plots to avoid confusion when changing methods.
 *
 * @example
 * setupSimulationMethodToggle(); // Activates input toggle and plot clearing
 */
function setupSimulationMethodToggle() {
  const simulationSelect = document.getElementById("simulation_method");
  const inputsMap = {
    NVT_Langevin_BAOAB: "inputs_nvt",
    muVT_Monte_Carlo: "inputs_muvt",
  };
  const plotIds = ["output", "positions3d", "temperature", "pressure", "plotMetaBias", "plotParticles3D", "plotProjection2D", "plotCVHistogram"];

  simulationSelect.addEventListener("change", () => {
    const method = simulationSelect.value;

    // Toggle inputs
    Object.entries(inputsMap).forEach(([key, id]) => {
      const el = document.getElementById(id);
      if (el) el.classList.toggle("hidden", key !== method);
    });

    // Limpiar gráficos
    plotIds.forEach((id) => {
      const gd = document.getElementById(id);
      if (!gd) return;
      Plotly.purge(gd);
      gd.innerHTML = "";
    });
  });
}


/**
 * Configura el toggle de modo oscuro en la interfaz.
 *
 * - Cambia la clase "dark" en el elemento <html> al hacer click en el botón.
 * - Actualiza el icono del botón (sol/luna) según el estado actual.
 * - Llama a `updatePlotlyTheme` para aplicar el tema oscuro a los gráficos.
 *
 * @example
 * setupDarkModeToggle(); // Activa el toggle de modo oscuro
 *
 * Sets up the dark mode toggle in the interface.
 *
 * - Toggles the "dark" class on the <html> element when the button is clicked.
 * - Updates the button icon (sun/moon) according to the current state.
 * - Calls `updatePlotlyTheme` to apply the dark theme to Plotly charts.
 *
 * @example
 * setupDarkModeToggle(); // Activates dark mode toggle
 */
function setupDarkModeToggle() {
  const toggleBtn = document.getElementById("dark-toggle");
  const icon = toggleBtn?.querySelector("i");

  if (!toggleBtn || !icon) return;

  const updateIcon = () => {
    const isDark = document.documentElement.classList.contains("dark");
    icon.classList.toggle("fa-moon", isDark);
    icon.classList.toggle("fa-sun", !isDark);
  };

  toggleBtn.addEventListener("click", () => {
    document.documentElement.classList.toggle("dark");
    updateIcon();
    updatePlotlyTheme();
  });

  updateIcon();
}


/**
 * Configura el modal de información con pestañas (tabs) y acordeones (accordions).
 *
 * Funcionalidades:
 * - Abrir y cerrar el modal.
 * - Cambiar entre pestañas y actualizar estilos de botones.
 * - Inicializar acordeones y alternar su estado.
 * - Renderizar MathJax dentro del modal si está disponible.
 *
 * @example
 * setupInfoModal(); // Inicializa el modal de información
 *
 * Sets up the info modal with tabs and accordions.
 *
 * Features:
 * - Open and close the modal.
 * - Switch between tabs and update button styles.
 * - Initialize accordions and toggle their state.
 * - Render MathJax inside the modal if available.
 *
 * @example
 * setupInfoModal(); // Initializes the info modal
 */
function setupInfoModal() {
  const infoIcon = document.querySelector("#infoModalTrigger");
  const infoModal = document.querySelector("#infoModal");
  const closeModal = document.querySelector("#closeModal");

  const tabButtons = [
    { button: document.getElementById("tab1Button"), content: document.getElementById("tab1Content") },
    { button: document.getElementById("tab2Button"), content: document.getElementById("tab2Content") },
  ];

  const setActiveTab = (activeIdx) => {
    tabButtons.forEach((tab, idx) => {
      const isActive = idx === activeIdx;
      tab.content.classList.toggle("hidden", !isActive);
      tab.button.classList.toggle("border-blue-600", isActive);
      tab.button.classList.toggle("dark:border-blue-400", isActive);
      tab.button.classList.toggle("text-blue-600", isActive);
      tab.button.classList.toggle("dark:text-blue-400", isActive);
      tab.button.classList.toggle("border-transparent", !isActive);
      tab.button.classList.toggle("text-gray-600", !isActive);
      tab.button.classList.toggle("dark:text-gray-300", !isActive);
    });
  };

  const resetAccordions = () => {
    tabButtons.forEach((tab) => {
      const accordions = tab.content.querySelectorAll(".accordion-content");
      const icons = tab.content.querySelectorAll(".accordion-btn i");
      accordions.forEach((content, idx) => content.classList.toggle("hidden", idx !== 0));
      icons.forEach((icon, idx) => icon.classList.toggle("rotate-180", idx === 0));
    });
  };

  const openModal = () => {
    infoModal.classList.remove("hidden");
    infoModal.classList.add("flex");
    resetAccordions();
    // Renderizar MathJax si existe
    window.MathJax?.typesetPromise([infoModal]);
  };

  const closeModalFn = () => {
    infoModal.classList.add("hidden");
    infoModal.classList.remove("flex");
    setActiveTab(0); // primera tab activa
    resetAccordions(); // acordeones al estado inicial
  };

  infoIcon?.addEventListener("click", openModal);
  closeModal?.addEventListener("click", closeModalFn);
  infoModal?.addEventListener("click", (e) => {
    if (e.target === infoModal) closeModalFn();
  });

  // Tabs
  tabButtons.forEach((tab, idx) => {
    tab.button.addEventListener("click", () => setActiveTab(idx));
  });

  // Acordeones
  document.querySelectorAll(".accordion-btn").forEach((btn) => {
    btn.addEventListener("click", () => {
      const content = btn.nextElementSibling;
      const icon = btn.querySelector("i");
      content.classList.toggle("hidden");
      icon.classList.toggle("rotate-180");
    });
  });
}


/**
 * Inicializa la interfaz de usuario (UI) de la aplicación.
 *
 * Funcionalidades incluidas:
 * 1. Carga el HTML del modal de información.
 * 2. Configura el overlay para dispositivos móviles o pantallas pequeñas.
 * 3. Configura el toggle de método de simulación y limpieza de gráficos.
 * 4. Configura el toggle de modo oscuro.
 * 5. Configura el modal de información con tabs y acordeones.
 *
 * @async
 * @returns {Promise<void>} - No devuelve ningún valor.
 *
 * @example
 * await initUI(); // Inicializa toda la interfaz de usuario
 *
 * Initializes the user interface (UI) of the application.
 *
 * Included functionalities:
 * 1. Loads the HTML for the info modal.
 * 2. Sets up the overlay for mobile devices or small screens.
 * 3. Sets up the simulation method toggle and plot clearing.
 * 4. Sets up the dark mode toggle.
 * 5. Sets up the info modal with tabs and accordions.
 *
 * @async
 * @returns {Promise<void>} - Does not return any value.
 *
 * @example
 * await initUI(); // Initializes the entire user interface
 */
export async function initUI() {
  await loadModalHTML();

  setupDeviceOverlay();
  setupSimulationMethodToggle();
  setupDarkModeToggle();
  setupInfoModal();
}
