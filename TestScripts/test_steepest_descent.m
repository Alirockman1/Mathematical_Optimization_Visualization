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
% objectiveFunction = @(x) (1 - x(:,1)).^2 + 100 * (x(:,2) - x(:,1).^2).^2;
% initialDesign = [0.1, 0.1];
% designSpaceLowerBound = [0,0];
% designSpaceUpperBound = [1.5, 1.5];

objectiveFunction = @(x) (x(:,1)-1.5).^2 + 5*(x(:,2)-2.0).^2;
initialDesign = [0.5, 1.5];
designSpaceLowerBound = [0,1];
designSpaceUpperBound = [3,3];

[optimumCandidate, optimumObjectiveValue, optimizationData] = optimization_steepest_descent(...
    objectiveFunction, initialDesign, designSpaceLowerBound, designSpaceUpperBound, ...
    'MaxIterations', 30, 'Tolerance', 1e-6, 'LineSearchFunction', @line_search_golden_ratio);

plot_point_optimization_results(optimizationData, 'SaveFolder',saveFolder,'SaveFigureOptions',{'Size',figureSize},...
    'LineWidth', 1, 'Color', 'b', 'Marker','o', 'LineStyle', '-');

plot_point_gradient_progress_slider(optimizationData,objectiveFunction, 'OptimizationMethod', 'Steepest descent', 'SaveFolder',saveFolder);


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;