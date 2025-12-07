# Green-Chemistry-Process-Simulator
A dynamic, first-principles kinetic simulator for the multi-metric analysis and visualisation of green chemistry principles in synthetic processes.

## Project Overview
This project implements a dynamic, first-principles kinetic simulation engine to evaluate chemical synthesis routes against the 12 Principles of Green Chemistry. Traditional green chemistry analysis relies on static, endpoint metrics calculated from completed reactions. This project develops a computational tool that simulates chemical processes from the ground up, providing a dynamic, time-resolved view of their 'greenness'.

By solving the underlying differential equations of reaction kinetics, this simulator goes beyond simply asking "How green is the final product?" to analyse "How does the efficiency, waste, and hazard profile of the reaction evolve over time?"

## Core Features
- **Generalised kinetic engine:** an object-oriented engine capable of simulating any user-defined chemical mechanism, including multi-step reactions with fast equilibria and reactive intermediates.
- **Dynamic metrics calculation:** a dedicated module that calculates a suite of key green chemistry metrics (e.g., Process Mass Intensity, Environmental Quotient) at every time step of the simulation.
- **Visualisation dashboard:** an interactive dashboard featuring an animated radar chart that visualises the evolving trade-offs of the process as it proceeds from reactants to products.
- **Rigorously validated:** the engine's accuracy is confirmed by a comprehensive validation suite, proving its adherence to fundamental physical laws including conservation of mass, thermodynamic equilibrium, and Arrhenius kinetics.

The tool is applied to a case study of aspirin synthesis, comparing traditional vs. greener routes through kinetic modelling rather than manual data transcription.

## Key Modules
### 1. The Kinetic Engine (`engine.py`)
A custom-built, object-oriented simulation engine in Python.
* **`ChemicalSpecies`:** encapsulates physical constants (including Van der Waals constants).
* **`Reaction`:** models elementary kinetic steps using the Arrhenius equation ($k = Ae^(-E_a/RT)$). Supports mass action kinetics.
* **`GeneralChemicalSystem`:** the core "solver" class. It constructs the system of ODEs ($dn/dt$) for an arbitrary reaction mechanism and solves them using the `scipy.integrate.solve_ivp` library (Radau method for stiff systems). It adheres strictly to the conservation of mass principle.

### 2. The Metrics Calculator (`metrics.py`)
A bridge between raw physical data (moles, temperature, time) and Green Chemistry insights.
* **Static metrics:** calculates standard industrial metrics including Atom Economy (AE) and static Process Mass Intensity (PMI).
* **Dynamic, time-series metrics:** uses vectorised NumPy operations to calculate the evolution of metrics over time:
  - *Dynamic PMI:* visualises how mass efficiency improves (or stalls) as conversion progresses.
  - *Dynamic hazard score:* monitors the reactor for transient spikes in toxicity (e.g., accumulation of hazardous intermediates) that endpoint analysis misses.
  - *Accumulated waste:* dissects the waste stream composition at every time step.

## Methodology & Novelty
Standard green chemistry metrics are often applied as a "post-mortem" label. This project introduces a paradigm shift to "live process diagnosis".
1. We treat "greenness" not as a fixed property, but as a dynamic variable. The simulator reveals when a process becomes inefficient (diminishing returns) or dangerous.
2. By building the model from fundamental kinetics ($A, E_a$) rather than empirical yield tables, the tool can predict the impact of changing conditions (temperature, concentration, catalyst loading) on the final environmental footprint.
3. The tool explicitly visualises the trade-offs between competing green principles (e.g., a route using safer reagents that suffers from poor atom economy or slow kinetics).
