function CalculateCommunicationProbabilitiesFromDeltaNotchOdeModel
%%% Script to explore how SoptSC calculates cell-cell communication
%%% probabilities for a juxtacrine signalling pathway, Delta-Notch. That
%%% is, Delta-Notch is a signalling pathway that is regulated by direct
%%% contact between adjacent cells. We generate synthetic Delta-Notch data
%%% by simulating the nonlinear ODE system described in Collier et al.
%%% (1996), using the parameter values from Germano and Osborne (2020), and
%%% then input these as gene expression values into SoptSC to see how
%%% SoptSC handles this sort of data.

clear all; close all; % Usual pre-amble

%%
% Set the data directory
dataDirectory = '~/Documents/scRNASeqAnalysisAndModelling/Julia/Solutions/';
dataFile = [dataDirectory, 'DeltaNotchOdeSolution.csv'];

% This solution is stored as a CSV, with headers t, N1, D1, .... Nm, Dm,
% where Ni and Di are the Notch and Delta expressions in cell i, where i
% ranges from 1 to m.
solution_data = readmatrix(dataFile); % Read in the CSV file
num_cells = 0.5 * (size(solution_data, 2) - 1); % Get the number of cells

% Get the Notch and Delta expressions at the beginning and end
initial_solution = solution_data(1, 2:end); % Initial solution; the first column is time, so we don't include that
initial_notch_expression = initial_solution(1:2:end); % Get the notch expression levels (alternates between notch and delta)
initial_delta_expression = initial_solution(2:2:end); % Get the delta expression lvels

final_solution = solution_data(end, 2:end); % The first column is time, so we don't include that
final_notch_expression = final_solution(1:2:end); % Get the notch expression levels (alternates between notch and delta)
final_delta_expression = final_solution(2:2:end); % Get the delta expression lvels

% We're plotting the solutions at the start and end time just to confirm
% they're doing what they're supposed to be doing. That is, when Delta's
% high, Notch should be low, and vice versa.
figure(1)
subplot(1, 2, 1)
hold on
plot(1:num_cells, initial_notch_expression)
plot(1:num_cells, initial_delta_expression)
title('Initial expression levels')

subplot(1, 2, 2)
hold on 
plot(1:num_cells, final_notch_expression)
plot(1:num_cells, final_delta_expression)
title('Final expression levels')

% Construct the gene expression data from the final solution
gene_expression_data = [final_notch_expression; final_delta_expression];
ligand = {'Delta'};
receptor = {'Notch'};
genes_in_data = {'Notch'; 'Delta'};

% Finally we can calculate the communication probabilities using SoptSC
[Pidv, Pall] = LR_Interaction(gene_expression_data, genes_in_data, ligand, receptor);

% Look at where the non-zero probabilities lie
figure(2)
spy(Pall)