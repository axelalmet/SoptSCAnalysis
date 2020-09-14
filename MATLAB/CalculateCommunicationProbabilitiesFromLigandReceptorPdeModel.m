function CalculateCommunicationProbabilitiesFromLigandReceptorPdeModel
%%% Script to explore how various situations explored in a 1D PDE model of
%%% ligand-receptor binding, as a proxy for cell-cell communication 
%%% We first look at the scenario when there's a ligand diffusing at one end only.
%%% We will assume that the bound receptor represents the downstream receptor
%%% expression within the cell. 
clear all; close all; % Clear all variables and close all figures

dataDirectory = '~/Documents/scRNASeqAnalysisAndModelling/PyPde/Results/LigandsSameRegion/';
%% First scenario: ligand at one end, constant free receptor throughout.
% First load the data.
dataFile = [dataDirectory, '1dLigandReceptorModelWithLigandProductionAtOneEndAndConstantReceptor.h5'];

% Structure is N x 3 x T, where N is the number of mesh points and T is the number of timesteps
pde_solutions = h5read(dataFile, '/data');

num_points = size(pde_solutions, 1); % Get the number of points from the mesh size
num_timepoints = size(pde_solutions, 3); % Get the number of timepoints

final_pde_solution = pde_solutions(:, :, num_timepoints); % Get the final solution at T = 5
mesh = linspace(0, 1, num_points); % Create the mesh from the number of points

figure(1)
subplot(3, 3, 1)
hold on
plot(mesh, final_pde_solution(:, 1)') % Plot the ligand
plot(mesh, final_pde_solution(:, 2)') % Plot the free receptor
plot(mesh, final_pde_solution(:, 3)') % Plot the bound ligand/receptor

% Let's try and calculate the cell-cell communication probabilities via
% SoptSC

ligand_receptor_data = final_pde_solution(:, [1 3])' + 1e-4; % Store the data
ligand = {'Wnt3'};
receptor = {'Fzd1'};
pde_genes = {'Wnt3'; 'Fzd1'};

% Finally we can calculate the communication probabilities using SoptSC
[Pidv, Pall] = LR_Interaction(ligand_receptor_data, pde_genes, ligand, receptor);

figure(2)
subplot(3, 3, 1)
spy(Pall)
%% Oscillating receptor level, defined such that every other cell has available receptor
% First load the data.
dataFile = [dataDirectory, '1dLigandReceptorModelWithLigandProductionAtOneEndAndOscillatingReceptor.h5'];

% Structure is N x 3 x T, where N is the number of mesh points and T is the number of timesteps
pde_solutions = h5read(dataFile, '/data'); 

num_points = size(pde_solutions, 1); 
num_timepoints = size(pde_solutions, 3);

final_pde_solution = pde_solutions(:, :, num_timepoints); % Get the final solution at T = 5
mesh = linspace(0, 1, num_points); % Create the mesh from the number of points

figure(1)
subplot(3, 3, 2)
hold on
plot(mesh, final_pde_solution(:, 1)') % Plot the ligand
plot(mesh, final_pde_solution(:, 2)') % Plot the free receptor
plot(mesh, final_pde_solution(:, 3)') % Plot the bound ligand/receptor

% Let's try and calculate the cell-cell communication probabilities via
% SoptSC

ligand_receptor_data = final_pde_solution(:, [1 3])' + 1e-4; % Add the small factor so that we don't get NaN during probability counts
ligand = {'Wnt3'};
receptor = {'Fzd1'};
pde_genes = {'Wnt3'; 'Fzd1'};

% Finally we can calculate the communication probabilities using SoptSC
[Pidv, Pall] = LR_Interaction(ligand_receptor_data, pde_genes, ligand, receptor);

figure(2)
subplot(3, 3, 2)
spy(Pall)
%% Slowly oscillating receptor level, so that only one every few cells expresses the required receptor
% First load the data.
dataFile = [dataDirectory, '1dLigandReceptorModelWithLigandProductionAtOneEndAndSlowlyOscillatingReceptor.h5'];

% Structure is N x 3 x T, where N is the number of mesh points and T is the number of timesteps
pde_solutions = h5read(dataFile, '/data'); 

num_points = size(pde_solutions, 1); 
num_timepoints = size(pde_solutions, 3);

final_pde_solution = pde_solutions(:, :, num_timepoints); % Get the final solution at T = 5
mesh = linspace(0, 1, num_points); % Create the mesh from the number of points

figure(1)
subplot(3, 3, 3)
hold on
plot(mesh, final_pde_solution(:, 1)') % Plot the ligand
plot(mesh, final_pde_solution(:, 2)') % Plot the free receptor
plot(mesh, final_pde_solution(:, 3)') % Plot the bound ligand/receptor

% Let's try and calculate the cell-cell communication probabilities via
% SoptSC

ligand_receptor_data = final_pde_solution(:, [1 3])' + 1e-4;  % Add the small factor so that we don't get NaN during probability counts
ligand = {'Wnt3'};
receptor = {'Fzd1'};
pde_genes = {'Wnt3'; 'Fzd1'};

% Finally we can calculate the communication probabilities using SoptSC
[Pidv, Pall] = LR_Interaction(ligand_receptor_data, pde_genes, ligand, receptor);

figure(2)
subplot(3, 3, 3)
spy(Pall)
%% Second scenario of the ligand from both ends: constant receptor availability
dataFile = [dataDirectory, '1dLigandReceptorModelWithLigandProductionAtBothEndsAndConstantReceptor.h5'];

% Structure is N x 3 x T, where N is the number of mesh points and T is the number of timesteps
pde_solutions = h5read(dataFile, '/data'); 

num_points = size(pde_solutions, 1); 
num_timepoints = size(pde_solutions, 3);

final_pde_solution = pde_solutions(:, :, num_timepoints); % Get the final solution at T = 5
mesh = linspace(0, 1, num_points); % Create the mesh from the number of points

figure(1)
subplot(3, 3, 4)
hold on
plot(mesh, final_pde_solution(:, 1)') % Plot the ligand
plot(mesh, final_pde_solution(:, 2)') % Plot the free receptor
plot(mesh, final_pde_solution(:, 3)') % Plot the bound ligand/receptor

% Let's try and calculate the cell-cell communication probabilities via
% SoptSC

ligand_receptor_data = final_pde_solution(:, [1 3])' + 1e-4; % Store the data
ligand = {'Wnt3'};
receptor = {'Fzd1'};
pde_genes = {'Wnt3'; 'Fzd1'};

% Finally we can calculate the communication probabilities using SoptSC
[Pidv, Pall] = LR_Interaction(ligand_receptor_data, pde_genes, ligand, receptor);

figure(2)
subplot(3, 3, 4)
spy(Pall)
%% Second scenario of the ligand from both ends: alternating receptor availability
dataFile = [dataDirectory, '1dLigandReceptorModelWithLigandProductionAtBothEndsAndOscillatingReceptor.h5'];

% Structure is N x 3 x T, where N is the number of mesh points and T is the number of timesteps
pde_solutions = h5read(dataFile, '/data'); 

num_points = size(pde_solutions, 1); 
num_timepoints = size(pde_solutions, 3);

final_pde_solution = pde_solutions(:, :, num_timepoints); % Get the final solution at T = 5
mesh = linspace(0, 1, num_points); % Create the mesh from the number of points

figure(1)
subplot(3, 3, 5)
hold on
plot(mesh, final_pde_solution(:, 1)') % Plot the ligand
plot(mesh, final_pde_solution(:, 2)') % Plot the free receptor
plot(mesh, final_pde_solution(:, 3)') % Plot the bound ligand/receptor

% Let's try and calculate the cell-cell communication probabilities via
% SoptSC

ligand_receptor_data = final_pde_solution(:, [1 3])' + 1e-4; % Store the data
ligand = {'Wnt3'};
receptor = {'Fzd1'};
pde_genes = {'Wnt3'; 'Fzd1'};

% Finally we can calculate the communication probabilities using SoptSC
[Pidv, Pall] = LR_Interaction(ligand_receptor_data, pde_genes, ligand, receptor);

figure(2)
subplot(3, 3, 5)
spy(Pall)
%% And again, ligand from both ends, but slowly oscillating receptor availability (i.e. more spatial distance)
dataFile = [dataDirectory, '1dLigandReceptorModelWithLigandProductionAtBothEndsAndSlowlyOscillatingReceptor.h5'];

% Structure is N x 3 x T, where N is the number of mesh points and T is the number of timesteps
pde_solutions = h5read(dataFile, '/data'); 

num_points = size(pde_solutions, 1); 
num_timepoints = size(pde_solutions, 3);

final_pde_solution = pde_solutions(:, :, num_timepoints); % Get the final solution at T = 5
mesh = linspace(0, 1, num_points); % Create the mesh from the number of points

figure(1)
subplot(3, 3, 6)
hold on
plot(mesh, final_pde_solution(:, 1)') % Plot the ligand
plot(mesh, final_pde_solution(:, 2)') % Plot the free receptor
plot(mesh, final_pde_solution(:, 3)') % Plot the bound ligand/receptor

% Let's try and calculate the cell-cell communication probabilities via
% SoptSC

ligand_receptor_data = final_pde_solution(:, [1 3])' + 1e-4; % Store the data
ligand = {'Wnt3'};
receptor = {'Fzd1'};
pde_genes = {'Wnt3'; 'Fzd1'};

% Finally we can calculate the communication probabilities using SoptSC
[Pidv, Pall] = LR_Interaction(ligand_receptor_data, pde_genes, ligand, receptor);

figure(2)
subplot(3, 3, 6)
spy(Pall)

%% Ligand production from middle, constant free receptor throughout.
% First load the data.
dataFile = [dataDirectory, '1dLigandReceptorModelWithLigandProductionInMiddleAndConstantReceptor.h5'];

% Structure is N x 3 x T, where N is the number of mesh points and T is the number of timesteps
pde_solutions = h5read(dataFile, '/data'); 

num_points = size(pde_solutions, 1); 
num_timepoints = size(pde_solutions, 3);

final_pde_solution = pde_solutions(:, :, num_timepoints); % Get the final solution at T = 5
mesh = linspace(0, 1, num_points); % Create the mesh from the number of points

figure(1)
subplot(3, 3, 7)
hold on
plot(mesh, final_pde_solution(:, 1)') % Plot the ligand
plot(mesh, final_pde_solution(:, 2)') % Plot the free receptor
plot(mesh, final_pde_solution(:, 3)') % Plot the bound ligand/receptor

% Let's try and calculate the cell-cell communication probabilities via
% SoptSC

ligand_receptor_data = final_pde_solution(:, [1 3])' + 1e-4; % Store the data
ligand = {'Wnt3'};
receptor = {'Fzd1'};
pde_genes = {'Wnt3'; 'Fzd1'};

% Finally we can calculate the communication probabilities using SoptSC
[Pidv, Pall] = LR_Interaction(ligand_receptor_data, pde_genes, ligand, receptor);

figure(2)
subplot(3, 3, 7)
spy(Pall)

%% Ligand production from middle, oscillating free receptor throughout.
% First load the data.
dataFile = [dataDirectory, '1dLigandReceptorModelWithLigandProductionInMiddleAndOscillatingReceptor.h5'];

% Structure is N x 3 x T, where N is the number of mesh points and T is the number of timesteps
pde_solutions = h5read(dataFile, '/data'); 

num_points = size(pde_solutions, 1); 
num_timepoints = size(pde_solutions, 3);

final_pde_solution = pde_solutions(:, :, num_timepoints); % Get the final solution at T = 5
mesh = linspace(0, 1, num_points); % Create the mesh from the number of points

figure(1)
subplot(3, 3, 8)
hold on
plot(mesh, final_pde_solution(:, 1)') % Plot the ligand
plot(mesh, final_pde_solution(:, 2)') % Plot the free receptor
plot(mesh, final_pde_solution(:, 3)') % Plot the bound ligand/receptor

% Let's try and calculate the cell-cell communication probabilities via
% SoptSC

ligand_receptor_data = final_pde_solution(:, [1 3])' + 1e-4; % Store the data
ligand = {'Wnt3'};
receptor = {'Fzd1'};
pde_genes = {'Wnt3'; 'Fzd1'};

% Finally we can calculate the communication probabilities using SoptSC
[Pidv, Pall] = LR_Interaction(ligand_receptor_data, pde_genes, ligand, receptor);

figure(2)
subplot(3, 3, 8)
spy(Pall)

%% Ligand production from middle, slowly oscillating free receptor throughout.
% First load the data.
dataFile = [dataDirectory, '1dLigandReceptorModelWithLigandProductionInMiddleAndSlowlyOscillatingReceptor.h5'];

% Structure is N x 3 x T, where N is the number of mesh points and T is the number of timesteps
pde_solutions = h5read(dataFile, '/data'); 

num_points = size(pde_solutions, 1); 
num_timepoints = size(pde_solutions, 3);

final_pde_solution = pde_solutions(:, :, num_timepoints); % Get the final solution at T = 5
mesh = linspace(0, 1, num_points); % Create the mesh from the number of points

figure(1)
subplot(3, 3, 9)
hold on
plot(mesh, final_pde_solution(:, 1)') % Plot the ligand
plot(mesh, final_pde_solution(:, 2)') % Plot the free receptor
plot(mesh, final_pde_solution(:, 3)') % Plot the bound ligand/receptor

% Let's try and calculate the cell-cell communication probabilities via
% SoptSC

ligand_receptor_data = final_pde_solution(:, [1 3])' + 1e-4; % Store the data
ligand = {'Wnt3'};
receptor = {'Fzd1'};
pde_genes = {'Wnt3'; 'Fzd1'};

% Finally we can calculate the communication probabilities using SoptSC
[Pidv, Pall] = LR_Interaction(ligand_receptor_data, pde_genes, ligand, receptor);

figure(2)
subplot(3, 3, 9)
spy(Pall)