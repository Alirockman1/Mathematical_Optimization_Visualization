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
objectiveFunction = @(x) -30 * x(:,1) -80 * x(:,2);
% inequalityConstraintFunction = @(x) [x(:,1) + x(:,2) - 18; x(:,1) + 2 * x(:,2) - 20; -x(:,1) + x(:,2) - 4];
inequalityConstraintFunction = {
    @(x) x(:,1) + x(:,2) - 18;
    @(x) x(:,1) + 2 * x(:,2) - 20;
    @(x) -x(:,1) + x(:,2) - 4};
equalityConstraintFunction = @(x) [];
initialDesign = [0.1, 0.1];
designSpaceLowerBound = [0,0];
designSpaceUpperBound = [10, 10];

[optimumCandidate, optimumObjectiveValue, optimizationData] = optimization_simplex(...
    objectiveFunction, inequalityConstraintFunction, equalityConstraintFunction, initialDesign, designSpaceLowerBound, designSpaceUpperBound, ...
    'MaxIterations', 30, 'Tolerance', 1e-6);

plot_point_optimization_results(optimizationData, 'SaveFolder',saveFolder,'SaveFigureOptions',{'Size',figureSize},...
     'LineWidth', 1, 'Color', 'b', 'Marker','o', 'LineStyle', '-');
 
plot_constraint_optimization_progress_slider(optimizationData, objectiveFunction, 'OptimizationMethod', 'Linear programming - Simplex', 'SaveFolder',saveFolder);

%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;