%% Cleanup
fclose all;
close all;
clear all;
clc;
more off;
diary off;


%% debugging
rng(4);


%% Documentation / Archive
rngState = rng;
saveFolder = save_diary_files(mfilename);
goldenRatio = (1+sqrt(5))/2;
figureSize = [goldenRatio 1]*8.5;


%% test
% sample multiObjective function
objectiveFunction = @(x) [(1 - x(:,1)).^2 , 100 * (x(:,2) - x(:,1).^2).^2];
% sample single objective function
% objectiveFunction = @(x) (1 - x(:,1)).^2 + 100 * (x(:,2) - x(:,1).^2).^2;
initialDesign = [0.1, 0.1];
designSpaceLowerBound = [0, 0];
designSpaceUpperBound = [1.5, 1.5];
constraints = @(x)deal(x(:,1)-0.5,[]);

options.BaseOptimizationFunction = @optimization_differential_evolution;
options.BaseOptimizationOptions = {
    'MaxIterations',100,...
    'PopulationSize',40,...
    'MutationFactor',0.85,...
    'CrossOverConstant',0.52
    };
options.ConvergenceCriterionFunction = @convergence_criterion_optimum_candidate_variance;
options.ConvergenceCriterionOptions = {};

[optimumCandidate, optimumObjectiveValue, optimizationData] = optimization_barrier_penalty(...
    objectiveFunction, initialDesign, designSpaceLowerBound, designSpaceUpperBound, constraints, options);

plot_point_optimization_results(optimizationData, 'SaveFolder',saveFolder,'SaveFigureOptions',{'Size',figureSize},...
    'LineWidth', 1, 'Color', 'b', 'Marker','o', 'LineStyle', '-');

% plot_point_multiObjective_optimization_results(optimizationData, 'SaveFolder',saveFolder,'SaveFigureOptions',{'Size',figureSize},...
%     'LineWidth', 1, 'Color', 'b', 'Marker','o', 'LineStyle', '-');

plot_contour_objective_constraint_results_slider(optimizationData, constraints, 'SaveFolder',saveFolder,'SaveFigureOptions',{'Size',figureSize});

%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

