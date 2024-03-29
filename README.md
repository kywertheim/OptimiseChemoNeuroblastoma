# OptimiseChemoNeuroblastoma
This repository contains two sets of code. One is used to simulate neuroblastoma's clonal evolution in the presence of vincristine and cyclophosphamide. The other is used to find the optimal chemotherapy schedule given an initial clonal composition.

Disclaimer: Matteo Italia wrote the simulation code. I coded the genetic algorithm, which Matteo modified.

Licence: GNU General Public License v3.0.

==========================================

Simulation.m.

Software requirement: MATLAB R2022a.

Execution: Running this m file within the MATLAB environment will initiate a deterministic simulation of neuroblastoma's clonal evolution in the presence of vincristine and cyclophosphamide.

Configuration: If the user wants to evaluate a different initial composition, they can change the initial percentage of cells that are resistant on line 57. If they want to change the distribution of resistant cells between the eight resistant clones, they can alter the code from line 126 to line 145. If they want to evaluate their own chemotherapy schedule, they can modify the code from line 82 to line 87.

Outputs: The program will generate a plot of the population dynamics automatically.

==========================================

GA.m.

Software requirement: MATLAB R2022a.

Execution: Running this m file within the MATLAB environment will find the optimal chemotherapy schedules for a broad range of initial clonal compositions.

Configuration: If the user wants to consider a different initial percentage of cells that are resistant, they can change the range on line 85. If they want to change the distribution of resistant cells between the eight resistant clones, they can alter the code from line 148 to line 167.

Outputs: The 12 solutions proposed for each initial clonal composition can be improved upon using the fmincon function, which searches the vicinity of each proposed solution.
