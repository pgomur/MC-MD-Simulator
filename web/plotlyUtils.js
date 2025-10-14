/**
 * Actualiza el tema de todos los gráficos de Plotly según el modo oscuro o claro.
 *
 * - Detecta si el documento tiene la clase "dark".
 * - Ajusta colores de fondo, ejes, títulos y leyenda de cada gráfico.
 * - Se aplica a todos los gráficos con IDs definidos en `plotIds`.
 *
 * @example
 * updatePlotlyTheme(); // Aplica el tema oscuro o claro a todos los gráficos existentes
 *
 * Updates the theme of all Plotly charts according to dark or light mode.
 *
 * - Detects if the document has the "dark" class.
 * - Adjusts background, axis, title, and legend colors for each chart.
 * - Applies to all charts with IDs defined in `plotIds`.
 *
 * @example
 * updatePlotlyTheme(); // Applies dark or light theme to all existing charts
 */
export function updatePlotlyTheme() {
  const isDark = document.documentElement.classList.contains("dark");
  const bgColor = isDark ? "#1f2937" : "#ffffff";
  const fontColor = isDark ? "#f9fafb" : "#111827";

  const plotIds = ["output", "positions3d", "temperature", "pressure", "plotMetaBias", "plotParticles3D", "plotProjection2D", "plotCVHistogram"];

  plotIds.forEach((id) => {
    const gd = document.getElementById(id);
    if (gd?.data) {
      Plotly.relayout(gd, {
        paper_bgcolor: bgColor,
        plot_bgcolor: bgColor,
        "xaxis.color": fontColor,
        "yaxis.color": fontColor,
        "title.font.color": fontColor,
        "legend.font.color": fontColor,
      });
    }
  });
}


/**
 * Genera un gráfico de energía (total, cinética y potencial) usando Plotly.
 *
 * - Convierte los pasos de simulación a tiempo (ps) y aplica factor de conversión de energía.
 * - Configura colores y fondo según el modo oscuro o claro.
 * - Muestra líneas con marcadores para cada tipo de energía.
 * - Configura leyenda, títulos de ejes y barras de herramientas.
 *
 * @param {Object} data - Datos de la simulación que contienen steps y energías.
 * @param {boolean} isDark - Indica si se debe aplicar el tema oscuro.
 *
 * @example
 * plotEnergy(simulationData, true); // Grafica la energía en modo oscuro
 *
 * Plots energy (total, kinetic, and potential) using Plotly.
 *
 * - Converts simulation steps to time (ps) and applies energy conversion factor.
 * - Sets colors and background based on dark or light mode.
 * - Shows lines with markers for each type of energy.
 * - Configures legend, axis titles, and toolbar.
 *
 * @param {Object} data - Simulation data containing steps and energy arrays.
 * @param {boolean} isDark - Indicates whether to apply dark mode theme.
 *
 * @example
 * plotEnergy(simulationData, false); // Plots energy in light mode
 */
export function plotEnergy(data, isDark) {
  const energyConversionFactor = 1;
  const timeConversionFactor = 2.15;
  const color = isDark ? "#e5e7eb" : "#000000";
  const bgColor = isDark ? "#1f2937" : "#ffffff";

  const energyTypes = [
    { key: "total_energy", name: "Total Energy" },
    { key: "kinetic_energy", name: "Kinetic Energy" },
    { key: "potential_energy", name: "Potential Energy" },
  ];

  const traces = energyTypes.map(({ key, name }) => ({
    x: data.steps.map((s) => s * timeConversionFactor),
    y: data[key].map((e) => e * energyConversionFactor),
    mode: "lines+markers",
    name,
    line: { width: 3 },
    marker: { size: 6 },
  }));

  Plotly.newPlot(
    "output",
    traces,
    {
      xaxis: { title: "Time (ps)", color },
      yaxis: { title: "Energía (kcal/mol)", color },
      margin: { l: 70, r: 70, t: 70, b: 70 },
      paper_bgcolor: bgColor,
      plot_bgcolor: bgColor,
      titlefont: { color },
      legend: {
        x: 0.5,
        y: 1.05,
        xanchor: "center",
        yanchor: "bottom",
        orientation: "h",
        font: { color },
      },
    },
    {
      responsive: true,
      displayModeBar: true,
      modeBarButtons: [["zoomIn2d", "zoomOut2d", "pan2d", "resetScale2d", "toImage"]],
    },
  );
}


/**
 * Genera un gráfico 3D animado de las posiciones de las partículas usando Plotly.
 *
 * Funcionalidades:
 * - Visualiza las posiciones X, Y, Z de cada partícula en cada paso de la simulación.
 * - Muestra la energía cinética por partícula en tooltips.
 * - Aplica colores a los marcadores según la coordenada Z.
 * - Incluye sliders para animar la simulación paso a paso.
 * - Se adapta a modo oscuro o claro.
 *
 * @param {Object} data - Datos de la simulación, incluyendo posiciones y energía cinética.
 * @param {number} n_particles - Número de partículas en la simulación.
 * @param {boolean} isDark - Indica si se debe aplicar el tema oscuro.
 *
 * @example
 * plotPositions3D(simulationData, 100, true); // Grafica posiciones 3D en modo oscuro
 *
 * Plots an animated 3D scatter of particle positions using Plotly.
 *
 * Features:
 * - Visualizes X, Y, Z positions of each particle at each simulation step.
 * - Shows kinetic energy per particle in hover tooltips.
 * - Colors markers based on Z-coordinate.
 * - Includes sliders to animate the simulation step by step.
 * - Adapts to dark or light mode.
 *
 * @param {Object} data - Simulation data including positions and kinetic energy.
 * @param {number} n_particles - Number of particles in the simulation.
 * @param {boolean} isDark - Indicates whether to apply dark mode theme.
 *
 * @example
 * plotPositions3D(simulationData, 50, false); // Plots 3D positions in light mode
 */
export function plotPositions3D(data, n_particles, isDark) {
  const color = isDark ? "#e5e7eb" : "#000000";
  const bgColor = isDark ? "#1f2937" : "#ffffff";
  const positions = Array.isArray(data.positions) ? data.positions : [];
  const kineticEnergy = Array.isArray(data.kinetic_energy) ? data.kinetic_energy : [];

  const createHoverText = (pos, kePerParticle) =>
    (pos?.[0] ?? []).map((x, idx) => {
      const y = pos[1]?.[idx]?.toFixed(2) ?? "N/A";
      const z = pos[2]?.[idx]?.toFixed(2) ?? "N/A";
      const keText = kePerParticle !== null ? `KE: ${kePerParticle.toFixed(3)}` : "KE: N/A";
      return `x: ${x.toFixed(2)}<br>y: ${y}<br>z: ${z}<br>${keText}`;
    });

  const createMarker = (zRaw = []) => {
    const zVals = zRaw.map(Number).filter(Number.isFinite);
    const hasZ = zVals.length > 0;
    return {
      size: 3,
      color: hasZ ? zRaw : "#FFD700",
      colorscale: [
        [0, "#00FFFF"],
        [0.25, "#00FF00"],
        [0.5, "#FFFF00"],
        [0.75, "#FFA500"],
        [1, "#FF4500"],
      ],
      cmin: hasZ ? Math.min(...zVals) : 0,
      cmax: hasZ ? Math.max(...zVals) : 1,
      opacity: 0.9,
    };
  };

  const frames = positions.map((pos, i) => {
    const kePerParticle = Number.isFinite(kineticEnergy[i]) && n_particles > 0 ? kineticEnergy[i] / n_particles : null;

    return {
      name: i.toString(),
      data: [
        {
          x: pos[0] ?? [],
          y: pos[1] ?? [],
          z: pos[2] ?? [],
          mode: "markers",
          type: "scatter3d",
          marker: createMarker(pos[2]),
          text: createHoverText(pos, kePerParticle),
          hoverinfo: "text",
        },
      ],
    };
  });

  const sliderSteps = frames.map((frame, i) => ({
    method: "animate",
    label: `Step ${i + 1}`,
    args: [[frame.name], { mode: "immediate", frame: { duration: 0, redraw: true }, transition: { duration: 0 } }],
  }));

  const layout3d = {
    margin: { l: 40, r: 40, b: 150, t: 0 },
    scene: {
      xaxis: { title: "X", color },
      yaxis: { title: "Y", color },
      zaxis: { title: "Z", color },
      aspectmode: "cube",
    },
    paper_bgcolor: bgColor,
    sliders: [
      {
        active: 0,
        currentvalue: { visible: true, font: { size: 14, color } },
        y: 0,
        steps: sliderSteps,
        tickcolor: "transparent",
        borderwidth: 0,
        bgcolor: isDark ? "#374151" : "#e5e7eb",
        font: { color: isDark ? "#f3f4f6" : "#111", size: 12 },
      },
    ],
    titlefont: { color },
  };

  const initialPos = positions[0] ?? [[], [], []];
  const initialKE = Number.isFinite(kineticEnergy[0]) && n_particles > 0 ? kineticEnergy[0] / n_particles : null;

  const trace3d = {
    x: initialPos[0],
    y: initialPos[1],
    z: initialPos[2],
    mode: "markers",
    type: "scatter3d",
    marker: createMarker(initialPos[2]),
    text: createHoverText(initialPos, initialKE),
    hoverinfo: "text",
  };

  Plotly.newPlot("positions3d", [trace3d], layout3d, {
    responsive: true,
    modeBarButtons: [["zoomIn2d", "zoomOut2d", "pan2d", "resetScale2d", "toImage"]],
  });

  if (frames.length) Plotly.addFrames("positions3d", frames);
}


/**
 * Genera un gráfico de temperatura de la simulación usando Plotly.
 *
 * Funcionalidades:
 * - Calcula la temperatura por paso a partir de la energía cinética y el número de partículas.
 * - Aplica factor de conversión de tiempo para el eje X.
 * - Configura colores y fondo según el modo oscuro o claro.
 * - Muestra líneas con marcadores para la temperatura.
 * - Incluye leyenda y herramientas de zoom/pan/exportación.
 *
 * @param {Object} data - Datos de la simulación, incluyendo steps y energía cinética.
 * @param {number} n_particles - Número de partículas en la simulación.
 * @param {boolean} isDark - Indica si se debe aplicar el tema oscuro.
 *
 * @example
 * plotTemperature(simulationData, 100, true); // Grafica la temperatura en modo oscuro
 *
 * Plots simulation temperature using Plotly.
 *
 * Features:
 * - Calculates temperature per step from kinetic energy and number of particles.
 * - Applies time conversion factor for the X-axis.
 * - Sets colors and background based on dark or light mode.
 * - Shows lines with markers for temperature.
 * - Includes legend and toolbar for zoom/pan/export.
 *
 * @param {Object} data - Simulation data including steps and kinetic energy.
 * @param {number} n_particles - Number of particles in the simulation.
 * @param {boolean} isDark - Indicates whether to apply dark mode theme.
 *
 * @example
 * plotTemperature(simulationData, 50, false); // Plots temperature in light mode
 */
export function plotTemperature(data, n_particles, isDark) {
  const timeConversionFactor = 2.15;
  const color = isDark ? "#e5e7eb" : "#000000";
  const bgColor = isDark ? "#1f2937" : "#ffffff";

  const temperatureY = (data.kinetic_energy ?? []).map((ke) => (n_particles > 0 && Number.isFinite(ke) ? (2 / (3 * n_particles)) * ke : null));

  const traceTemperature = {
    x: (data.steps ?? []).map((s) => s * timeConversionFactor),
    y: temperatureY,
    mode: "lines+markers",
    name: "Temperature",
    line: { width: 3 },
    marker: { size: 6 },
  };

  const layout = {
    xaxis: { title: "Time (ps)", color },
    yaxis: { title: "Temperature (ε/k_B)", color },
    margin: { l: 70, r: 70, t: 70, b: 70 },
    paper_bgcolor: bgColor,
    plot_bgcolor: bgColor,
    titlefont: { color },
    showlegend: true,
    legend: {
      x: 0.5,
      y: 1.05,
      xanchor: "center",
      yanchor: "bottom",
      orientation: "h",
      font: { color },
    },
  };

  const config = {
    responsive: true,
    displayModeBar: true,
    modeBarButtons: [["zoomIn2d", "zoomOut2d", "pan2d", "resetScale2d", "toImage"]],
  };

  Plotly.newPlot("temperature", [traceTemperature], layout, config);
}


/**
 * Genera un gráfico de presión de la simulación usando Plotly.
 *
 * Funcionalidades:
 * - Calcula la presión por paso usando la ecuación de estado: P = (N*T)/V + (r·F)/(3V)
 *   donde r·F es el producto escalar de posiciones y fuerzas.
 * - Convierte los pasos de simulación a tiempo (ps) para el eje X.
 * - Aplica colores y fondo según el modo oscuro o claro.
 * - Muestra líneas con marcadores para la presión.
 * - Incluye leyenda y herramientas de zoom/pan/exportación.
 *
 * @param {Object} data - Datos de la simulación, incluyendo posiciones, fuerzas, energía cinética y steps.
 * @param {number} n_particles - Número de partículas en la simulación.
 * @param {number} box_size - Tamaño de la caja cúbica de la simulación.
 * @param {boolean} isDark - Indica si se debe aplicar el tema oscuro.
 *
 * @example
 * plotPressure(simulationData, 100, 10, true); // Grafica la presión en modo oscuro
 *
 * Plots simulation pressure using Plotly.
 *
 * Features:
 * - Calculates pressure per step using the equation: P = (N*T)/V + (r·F)/(3V)
 *   where r·F is the dot product of positions and forces.
 * - Converts simulation steps to time (ps) for the X-axis.
 * - Sets colors and background based on dark or light mode.
 * - Shows lines with markers for pressure.
 * - Includes legend and toolbar for zoom/pan/export.
 *
 * @param {Object} data - Simulation data including positions, forces, kinetic energy, and steps.
 * @param {number} n_particles - Number of particles in the simulation.
 * @param {number} box_size - Cubic box size of the simulation.
 * @param {boolean} isDark - Indicates whether to apply dark mode theme.
 *
 * @example
 * plotPressure(simulationData, 50, 8, false); // Plots pressure in light mode
 */
export function plotPressure(data, n_particles, box_size, isDark) {
  const timeConversionFactor = 2.15;
  const color = isDark ? "#e5e7eb" : "#000000";
  const bgColor = isDark ? "#1f2937" : "#ffffff";
  const V = box_size > 0 ? Math.pow(box_size, 3) : null;

  const pressureData = (data.positions ?? []).map((pos, frameIdx) => {
    if (!pos?.[0] || !V || n_particles <= 0) return null;

    const keFrame = data.kinetic_energy?.[frameIdx];
    const T = Number.isFinite(keFrame) ? (2 / (3 * n_particles)) * keFrame : 0;

    const forcesFrame = data.forces?.[frameIdx];
    let rDotF = 0;

    if (forcesFrame) {
      let fxArr, fyArr, fzArr;

      if (forcesFrame.length >= 3 && forcesFrame[0].length === n_particles) {
        [fxArr, fyArr, fzArr] = forcesFrame;
      } else if (forcesFrame.length === n_particles && forcesFrame[0].length >= 3) {
        fxArr = forcesFrame.map((f) => f[0]);
        fyArr = forcesFrame.map((f) => f[1]);
        fzArr = forcesFrame.map((f) => f[2]);
      }

      if (fxArr && fyArr && fzArr) {
        for (let i = 0; i < n_particles; i++) {
          rDotF += (pos[0][i] ?? 0) * (fxArr[i] ?? 0) + (pos[1][i] ?? 0) * (fyArr[i] ?? 0) + (pos[2][i] ?? 0) * (fzArr[i] ?? 0);
        }
      }
    }

    const pressure = (n_particles * T) / V + rDotF / (3 * V);
    return Number.isFinite(pressure) ? pressure : null;
  });

  const tracePressure = {
    x: (data.steps ?? []).map((s) => s * timeConversionFactor),
    y: pressureData,
    mode: "lines+markers",
    name: "Pressure",
    line: { width: 3 },
    marker: { size: 6 },
  };

  const layout = {
    xaxis: { title: "Time (ps)", color },
    yaxis: { title: "Pressure (ε/σ³)", color },
    margin: { l: 70, r: 70, t: 70, b: 70 },
    paper_bgcolor: bgColor,
    plot_bgcolor: bgColor,
    titlefont: { color },
    showlegend: true,
    legend: {
      x: 0.5,
      y: 1.05,
      xanchor: "center",
      yanchor: "bottom",
      orientation: "h",
      font: { color },
    },
  };

  const config = {
    responsive: true,
    displayModeBar: true,
    modeBarButtons: [["zoomIn2d", "zoomOut2d", "pan2d", "resetScale2d", "toImage"]],
  };

  Plotly.newPlot("pressure", [tracePressure], layout, config);
}

// -----------------------------------------------------------------------------------------------------------------------------

/**
 * Genera un gráfico de Meta Bias (sesgo metadínamico) usando Plotly.
 *
 * Funcionalidades:
 * - Grafica los valores de meta_bias en función del índice del bin de la CV.
 * - Aplica colores y fondo según el modo oscuro o claro.
 * - Muestra líneas con marcadores y hover con formato personalizado.
 * - Incluye leyenda y herramientas de zoom/pan/exportación.
 *
 * @param {Object} data - Datos de la simulación, incluyendo el array meta_bias.
 * @param {boolean} isDark - Indica si se debe aplicar el tema oscuro.
 *
 * @example
 * plotMetaBias(simulationData, true); // Grafica Meta Bias en modo oscuro
 *
 * Plots Meta Bias (metadynamics bias) using Plotly.
 *
 * Features:
 * - Plots meta_bias values against the CV bin index.
 * - Applies colors and background based on dark or light mode.
 * - Shows lines with markers and custom hover text.
 * - Includes legend and toolbar for zoom/pan/export.
 *
 * @param {Object} data - Simulation data including meta_bias array.
 * @param {boolean} isDark - Indicates whether to apply dark mode theme.
 *
 * @example
 * plotMetaBias(simulationData, false); // Plots Meta Bias in light mode
 */
export function plotMetaBias(data, isDark) {
  const color = isDark ? "#e5e7eb" : "#000000";
  const bgColor = isDark ? "#1f2937" : "#ffffff";

  const trace = {
    x: (data.meta_bias ?? []).map((_, i) => i),
    y: data.meta_bias ?? [],
    mode: "lines+markers",
    type: "scatter",
    name: "Meta Bias",
    line: { color: isDark ? "#4ade80" : "#166534", width: 3 },
    marker: { size: 6, color: isDark ? "#22c55e" : "#14532d" },
    hovertemplate: "Bin CV: %{x} (Dimensionless)<br>Bias: %{y:.3f} kBT<extra></extra>",
  };

  const layout = {
    xaxis: { title: "CV bin index (dimensionless)", color },
    yaxis: { title: "Bias (kBT)", color },
    margin: { l: 70, r: 70, t: 70, b: 70 },
    paper_bgcolor: bgColor,
    plot_bgcolor: bgColor,
    titlefont: { color },
    showlegend: true,
    legend: {
      x: 0.5,
      y: 1.05,
      xanchor: "center",
      yanchor: "bottom",
      orientation: "h",
      font: { color },
    },
  };

  const config = {
    responsive: true,
    displayModeBar: true,
    modeBarButtons: [["zoomIn2d", "zoomOut2d", "pan2d", "resetScale2d", "toImage"]],
  };

  Plotly.newPlot("plotMetaBias", [trace], layout, config);
}


/**
 * Genera un gráfico 3D de posiciones de partículas usando Plotly.
 *
 * Funcionalidades:
 * - Representa las posiciones X, Y, Z de las partículas en 3D.
 * - Aplica colores a las partículas según su coordenada Z o un color fijo si no hay Z.
 * - Ajusta los colores, fuente y fondo según el modo oscuro o claro.
 * - Muestra hover con coordenadas de cada partícula.
 * - Incluye herramientas de zoom, paneo, reinicio de cámara y exportación.
 *
 * @param {Object} data - Datos de la simulación, incluyendo posiciones: [x[], y[], z[]].
 * @param {boolean} isDark - Indica si se debe aplicar el tema oscuro.
 *
 * @example
 * plotParticles3D(simulationData, true); // Grafica posiciones en modo oscuro
 *
 * Plots a 3D scatter of particle positions using Plotly.
 *
 * Features:
 * - Displays X, Y, Z positions of particles in 3D.
 * - Colors particles by Z coordinate or fixed color if Z is not available.
 * - Adjusts colors, font, and background according to dark or light mode.
 * - Shows hover info with particle coordinates.
 * - Includes zoom, pan, reset camera, and export tools.
 *
 * @param {Object} data - Simulation data, including positions: [x[], y[], z[]].
 * @param {boolean} isDark - Indicates whether to apply dark mode theme.
 *
 * @example
 * plotParticles3D(simulationData, false); // Plots 3D particle positions in light mode
 */
export function plotParticles3D(data, isDark) {
  const colorFont = isDark ? "#f9fafb" : "#111827";
  const bgColor = isDark ? "#1f2937" : "#ffffff";

  const titleElement = document.getElementById("plotTitle");
  if (titleElement && data.n_particles !== undefined) {
    titleElement.textContent = `3D Particle Distribution (N = ${data.n_particles})`;
  }

  const x = data.positions?.[0] ?? [];
  const y = data.positions?.[1] ?? [];
  const zRaw = data.positions?.[2] ?? [];

  const zVals = zRaw.map(Number).filter(Number.isFinite);
  const hasZ = zVals.length > 0;
  const zMin = hasZ ? Math.min(...zVals) : 0;
  const zMax = hasZ ? Math.max(...zVals) : 1;

  const marker = {
    size: 3,
    color: hasZ ? zRaw : "#FFD700",
    colorscale: [
      [0, "#00FFFF"],
      [0.25, "#00FF00"],
      [0.5, "#FFFF00"],
      [0.75, "#FFA500"],
      [1, "#FF4500"],
    ],
    cmin: zMin,
    cmax: zMax,
    opacity: 0.9,
  };

  const trace = {
    x,
    y,
    z: zRaw,
    mode: "markers",
    type: "scatter3d",
    marker,
    hovertemplate: "X: %{x:.2f} Å<br>Y: %{y:.2f} Å<br>Z: %{z:.2f} Å<extra></extra>",
  };

  const layout = {
    scene: {
      xaxis: { title: "X (Å)", color: colorFont },
      yaxis: { title: "Y (Å)", color: colorFont },
      zaxis: { title: "Z (Å)", color: colorFont },
      aspectmode: "cube",
    },
    paper_bgcolor: bgColor,
    plot_bgcolor: bgColor,
    font: { color: colorFont },
    margin: { l: 0, r: 0, b: 0, t: 0 },
  };

  Plotly.newPlot("plotParticles3D", [trace], layout, {
    responsive: true,
    modeBarButtons: [["zoom3d", "pan3d", "resetCameraDefault3d", "toImage"]],
  });
}


/**
 * Genera un gráfico 2D de proyección de posiciones de partículas usando Plotly.
 *
 * Funcionalidades:
 * - Proyecta las posiciones de las partículas en el plano XY.
 * - Aplica colores y fondo según el modo oscuro o claro.
 * - Muestra hover con coordenadas X e Y de cada partícula.
 * - Incluye herramientas de zoom, paneo y exportación.
 *
 * @param {Object} data - Datos de la simulación, incluyendo posiciones: [x[], y[], z[]].
 * @param {boolean} isDark - Indica si se debe aplicar el tema oscuro.
 *
 * @example
 * plotProjection2D(simulationData, true); // Grafica proyección XY en modo oscuro
 *
 * Plots a 2D projection of particle positions using Plotly.
 *
 * Features:
 * - Projects particle positions onto the XY plane.
 * - Adjusts colors and background based on dark or light mode.
 * - Shows hover info with X and Y coordinates.
 * - Includes zoom, pan, and export tools.
 *
 * @param {Object} data - Simulation data, including positions: [x[], y[], z[]].
 * @param {boolean} isDark - Indicates whether to apply dark mode theme.
 *
 * @example
 * plotProjection2D(simulationData, false); // Plots XY projection in light mode
 */
export function plotProjection2D(data, isDark) {
  const colorFont = isDark ? "#f9fafb" : "#111827";
  const bgColor = isDark ? "#1f2937" : "#ffffff";
  const markerColor = isDark ? "#f472b6" : "#9d174d";

  const trace = {
    x: data.positions?.[0] ?? [],
    y: data.positions?.[1] ?? [],
    mode: "markers",
    type: "scatter",
    marker: { size: 4, color: markerColor },
    hovertemplate: "X: %{x:.2f} Å<br>Y: %{y:.2f} Å<extra></extra>",
  };

  const layout = {
    xaxis: { title: "X (Å)", color: colorFont },
    yaxis: { title: "Y (Å)", color: colorFont },
    paper_bgcolor: bgColor,
    plot_bgcolor: bgColor,
    font: { color: colorFont },
    margin: { l: 70, r: 70, t: 50, b: 70 },
    titlefont: { color: colorFont },
    showlegend: false,
    legend: {
      x: 0.5,
      y: 1.05,
      xanchor: "center",
      yanchor: "bottom",
      orientation: "h",
      font: { color: colorFont },
    },
  };

  Plotly.newPlot("plotProjection2D", [trace], layout, {
    responsive: true,
    modeBarButtons: [["zoom3d", "pan3d", "resetCameraDefault3d", "toImage"]],
  });
}


/**
 * Genera un histograma del Bias de la CV (meta-bias) usando Plotly.
 *
 * Funcionalidades:
 * - Grafica la distribución de valores de meta_bias.
 * - Aplica colores y fondo según el modo oscuro o claro.
 * - Muestra hover con valor del bias y frecuencia.
 * - Incluye herramientas de zoom, paneo y exportación.
 *
 * @param {Object} data - Datos de la simulación, incluyendo meta_bias (array de sesgos).
 * @param {boolean} isDark - Indica si se debe aplicar el tema oscuro.
 *
 * @example
 * plotCVHistogram(simulationData, true); // Grafica histograma en modo oscuro
 *
 * Plots a histogram of CV bias (meta-bias) using Plotly.
 *
 * Features:
 * - Shows distribution of meta_bias values.
 * - Adjusts colors and background based on dark or light mode.
 * - Shows hover info with bias value and frequency.
 * - Includes zoom, pan, and export tools.
 *
 * @param {Object} data - Simulation data including meta_bias (array of bias values).
 * @param {boolean} isDark - Indicates whether to apply dark mode theme.
 *
 * @example
 * plotCVHistogram(simulationData, false); // Plots histogram in light mode
 */
export function plotCVHistogram(data, isDark) {
  const colorFont = isDark ? "#f9fafb" : "#111827";
  const bgColor = isDark ? "#1f2937" : "#ffffff";
  const markerColor = isDark ? "#facc15" : "#92400e";

  const trace = {
    x: data.meta_bias ?? [],
    type: "histogram",
    marker: { color: markerColor },
    hovertemplate: "Bias: %{x:.3f} kBT<br>Frequency: %{y}<extra></extra>",
  };

  const layout = {
    xaxis: { title: "Bias (kBT)", color: colorFont },
    yaxis: { title: "Frequency (number of bins)", color: colorFont },
    paper_bgcolor: bgColor,
    plot_bgcolor: bgColor,
    font: { color: colorFont },
    margin: { l: 70, r: 70, t: 50, b: 70 },
    titlefont: { color: colorFont },
    showlegend: false,
    legend: {
      x: 0.5,
      y: 1.05,
      xanchor: "center",
      yanchor: "bottom",
      orientation: "h",
      font: { color: colorFont },
    },
  };

  Plotly.newPlot("plotCVHistogram", [trace], layout, {
    responsive: true,
    modeBarButtons: [["zoom3d", "pan3d", "resetCameraDefault3d", "toImage"]],
  });
}
