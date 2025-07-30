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
% sample single objective function
objectiveFunction = @(x) (1 - x(:,1)).^2 + 100 * (x(:,2) - x(:,1).^2).^2;
% sample multiObjective function
% objectiveFunction = @(x) [(1 - x(:,1)).^2 , 100 * (x(:,2) - x(:,1).^2).^2];
initialDesign = [0.1, 0.1];
designSpaceLowerBound = [0, 0];
designSpaceUpperBound = [1.5, 1.5];

[optimumCandidate, optimumObjectiveValue, optimizationData] = optimization_differential_evolution(...
    objectiveFunction, initialDesign, designSpaceLowerBound, designSpaceUpperBound,...
    'MaxIterations', 100, 'PopulationSize', 10, 'MutationFactor', 0.85, 'CrossOverConstant', 0.52);

plot_point_optimization_results(optimizationData, 'SaveFolder',saveFolder,'SaveFigureOptions',{'Size',figureSize},...
    'LineWidth', 1, 'Color', 'b', 'Marker','o', 'LineStyle', '-');

% plot_point_multiObjective_optimization_results(optimizationData, 'SaveFolder',saveFolder,'SaveFigureOptions',{'Size',figureSize},...
%     'LineWidth', 1, 'Color', 'b', 'Marker','o', 'LineStyle', '-');
 
plot_differential_evolution_progress_nd_slider(optimizationData,'SaveFolder',saveFolder);

if length(optimizationData.IterationData(1).OptimumCandidate) == 2
    plot_objective_value_progress_2d(optimizationData,'Algorithm','Differential evolution','SaveFolder',saveFolder);
end

%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

