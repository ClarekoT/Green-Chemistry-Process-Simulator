# Green-Chemistry-Process-Simulator
A dynamic, first-principles kinetic simulator for the multi-metric analysis and visualisation of green chemistry principles in synthetic processes.

## Project Overview
Traditional green chemistry analysis relies on static, endpoint metrics calculated from completed reactions. This project develops a computational tool that simulates chemical processes from the ground up, providing a dynamic, time-resolved view of their 'greenness'.

By solving the underlying differential equations of reaction kinetics, this simulator goes beyond simply asking "How green is the final product?" to analyse "How does the efficiency, waste, and hazard profile of the reaction evolve over time?"

## Core Features
- **Generalised kinetic engine:** an object-oriented engine capable of simulating any user-defined chemical mechanism, including multi-step reactions with fast equilibria and reactive intermediates.
- **Dynamic metrics calculation:** a dedicated module that calculates a suite of key green chemistry metrics (e.g., Process Mass Intensity, Environmental Quotient) at every time step of the simulation.
- **Visualisation dashboard:** an interactive dashboard featuring an animated radar chart that visualises the evolving trade-offs of the process as it proceeds from reactants to products.
- **Rigorously validated:** the engine's accuracy is confirmed by a comprehensive validation suite, proving its adherence to fundamental physical laws including conservation of mass, thermodynamic equilibrium, and Arrhenius kinetics.

The tool is applied to a case study of aspirin synthesis, comparing traditional vs. greener routes through kinetic modelling rather than manual data transcription.
