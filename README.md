# 🧠 MATLAB Genetic Algorithm for Hybrid Energy System Optimization

This repository provides a **MATLAB-based simulation and optimization tool** that uses a **Genetic Algorithm (GA)** to determine the optimal configuration of a hybrid energy system. The system includes solar panels, wind turbines, batteries, and diesel generators, designed to meet dynamic hourly energy demand while minimizing the **Levelized Cost of Energy (LCOE)**.

---

## 🚀 Features

- **Optimization Algorithm:**
  - Custom-built Genetic Algorithm for multi-component energy systems.
  - Includes selection, crossover, and mutation operators.
  - Convergence criteria based on LCOE error threshold and generation count.

- **Energy System Modeling:**
  - Supports solar PV, wind turbines, diesel generators, and battery banks.
  - Realistic power output calculations based on weather data.
  - SOC (State of Charge) tracking for battery modeling.

- **Smart Resource Allocation:**
  - Adapts to constraints like available installation area.
  - Uses probabilistic mutation and selection for diverse solution exploration.

- **Interactive Progress Visualization:**
  - Dynamic progress bar with estimated time remaining.
  - Generates post-optimization plots and pie charts summarizing energy contributions.

- **Modular Structure:**
  - Key logic separated into functions:
    - `algoritmoGenetico.m`: Main optimization controller.
    - `planta_new.m`: Energy calculation and simulation.
    - `reproduccion.m`: GA crossover logic.
    - `mutacion.m`: GA mutation logic.
    - `specialTurbine.m`: Wind turbine performance modeling.
    - `cargaExcel.m`: Input parameter loading.