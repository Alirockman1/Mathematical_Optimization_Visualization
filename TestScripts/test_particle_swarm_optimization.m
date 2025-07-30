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
objectiveFunction = @(x) (1 - x(:,1)).^2 + 100 * (x(:,2) - x(:,1).^2).^2;
%objectiveFunction = @(x) [(1 - x(:,1)).^2 , 100 * (x(:,2) - x(:,1).^2).^2];
initialDesign = [0.1, 0.1];
designSpaceLowerBound = [0,0];
designSpaceUpperBound = [1.5, 1.5];
constraints = @(x)deal(x(:,1)-0.5,[]);

%% test - 4 dimensions
%objectiveFunction = @(x) (1 - x(:,1)).^2 + 100 * (x(:,2) - x(:,1).^2).^2;
% objectiveFunction = @(x) sum(100*(x(:,2:end)-x(:,1:end-1).^2).^2 + (1-x(:,1:end-1)).^2, 2);
% initialDesign = [0.1, 0.1, 0.1, 0.1];
% designSpaceLowerBound = [0, 0, 0, 0];
% designSpaceUpperBound = [1.5, 1.5, 1.5, 1.5];

[optimumCandidate, optimumObjectiveValue, optimizationData] = optimization_particle_swarm(...
    objectiveFunction,initialDesign,designSpaceLowerBound, designSpaceUpperBound,...
    'MaxIterations',100, 'PopulationSize', 10, ...
    'InertiaCoefficient',0.95,'CognitiveCoefficient',0.9,'SocialCoefficient',1.9);

plot_point_optimization_results(optimizationData, 'SaveFolder',saveFolder,'SaveFigureOptions',{'Size',figureSize},...
    'LineWidth', 1, 'Color', 'b', 'Marker','o', 'LineStyle', '-');

% plot_point_multiObjective_optimization_results(optimizationData, 'SaveFolder',saveFolder,'SaveFigureOptions',{'Size',figureSize},...
%     'LineWidth', 1, 'Color', 'b', 'Marker','o', 'LineStyle', '-');

plot_particle_swarm_progress_nd_slider(optimizationData,'SaveFolder',saveFolder);

if length(optimizationData.IterationData(1).OptimumCandidate) == 2
    plot_objective_value_progress_2d(optimizationData,'Algorithm','Particle Swarm','SaveFolder',saveFolder);
end


%% Save and Stop Transcripting
save([saveFolder,'Data.mat']);
diary off;

