function plot_point_multiObjective_optimization_results(optimizationData, varargin)
% PLOT_POINT_MULTIOBJECTIVE_OPTIMIZATION_RESULTS Visualizes results from multi-objective optimization
%
%   PLOT_POINT_MULTIOBJECTIVE_OPTIMIZATION_RESULTS(OPTIMIZATIONDATA) creates 
%   multiple visualization plots for multi-objective optimization results.
%
%   PLOT_POINT_MULTIOBJECTIVE_OPTIMIZATION_RESULTS(..., 'PARAM', VALUE) 
%   specifies additional parameters for customizing the plots.
%
%   INPUTS:
%       optimizationData - Structure containing optimization results with fields:
%           .ProblemData - Problem definition including bounds and options
%           .IterationData - Array of structs with optimization history
%           .InitialData - Initial conditions and values
%
%   OPTIONAL PARAMETERS:
%       'SaveFolder' - Path to folder for saving figures (default: [])
%       'CloseFigureAfterSaving' - Close figures after saving (default: false)
%       'SaveFigureOptions' - Cell array of options for figure saving
%       'LineWidth' - Line width for plots (default: 2)
%       'Color' - Line color (default: 'b')
%       'MarkerEdgeColor' - Marker edge color (default: 'b')
%       'Marker' - Marker style (default: 'o')
%       'LineStyle' - Line style (default: '-')
%       'LineType' - Type of line plot (default: 'line')
%       'MarkerFaceColor' - Marker fill color (default: 'auto')
%       'MarkerSize' - Marker size (default: 36)
%       'AxisFontSize' - Font size for axis labels (default: 12)
%       'ShowGrid' - Show grid lines ('on' or 'off', default: 'on')
%       'DesignVariableNames' - Cell array of names for design variables
%       'ObjectiveValueNames' - Cell array of names for objectives
%       'UseDesignSpaceAsAxisLimits' - Use design bounds for axis limits (default: false)
%       'PlotContour2DUseDesignSpaceLimit' - Use bounds for contour plots (default: false)
%       'PlotContour2DOptions' - Additional options for contour plots
%
%   OUTPUTS:
%       Generates multiple figures showing optimization progress:
%       1. Pairwise scatter matrix of objective values
%       2. Individual objective convergence plots
%       3. Normalized objective convergence
%       4. Error convergence plots
%
%   DETAILED PLOTS:
%       - Pairwise scatter matrix showing relationships between objectives
%       - Individual convergence plots for each objective function
%       - Normalized versions showing relative progress
%       - Error plots showing distance to final solution
%
%   SEE ALSO:
%       optimization_sequential_quadratic_programming, save_print_figure
%
%   Copyright 2025 Eduardo Rodrigues Della Noce
%   SPDX-License-Identifier: Apache-2.0

%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
% 
%       http://www.apache.org/licenses/LICENSE-2.0
% 
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

    p = inputParser;
    addParameter(p, 'SaveFolder',[]);
    addParameter(p, 'CloseFigureAfterSaving',false);
    addParameter(p, 'SaveFigureOptions',{});
    addParameter(p, 'LineWidth', 2, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'Color', 'b', @(x) ischar(x) || isstring(x));
    addParameter(p, 'MarkerEdgeColor', 'b', @(x) ischar(x) || isstring(x));
    addParameter(p, 'Marker', 'o', @(x) ischar(x) || isstring(x));
    addParameter(p, 'LineStyle', '-', @(x) ischar(x) || isstring(x));
    addParameter(p, 'LineType', 'line', @(x) ischar(x) || isstring(x));
    addParameter(p, 'MarkerFaceColor', 'auto', @(x) ischar(x) || isstring(x));
    addParameter(p, 'MarkerSize', 36, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'AxisFontSize', 12, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'ShowGrid', 'on', @(x) islogical(x) || ischar(x) && any(strcmpi(x, {'on', 'off'})));
    addParameter(p, 'DesignVariableNames',{});
    addParameter(p, 'ObjectiveValueNames',{});
    addParameter(p, 'UseDesignSpaceAsAxisLimits',false);
    addParameter(p, 'PlotContour2DUseDesignSpaceLimit',false);
    addParameter(p, 'PlotContour2DOptions',{});
    
    parse(p, varargin{:});
    options = p.Results;
    
    % get problem information / final solution
    nDimension = size(optimizationData.ProblemData.DesignSpaceLowerBound,2);
    optimalDesign = optimizationData.IterationData(end).OptimumCandidate;
    optimalObjectiveValue = optimizationData.IterationData(end).OptimumCandidateObjectiveValue;

    % extract optimum data
    optimumCandidate = [optimizationData.InitialData.OptimumCandidate;...
        vertcat(optimizationData.IterationData(:).OptimumCandidate)];
    optimumCandidateObjectiveValue = [optimizationData.InitialData.OptimumCandidateObjectiveValue;...
        vertcat(optimizationData.IterationData(:).OptimumCandidateObjectiveValue)];
    iterationIndex = 0:length(optimizationData.IterationData);

    % storing multiObjective values depending on unconstrained/constrained objective value
    if isfield(optimizationData.ProblemData.Options, 'BaseOptimizationFunction')
        for i=1:length(optimizationData.IterationData)
            optimumCandidateBaseObjectiveValue(i,:) = optimizationData.IterationData(i).BaseOptimizationData.IterationData(end).OptimumCandidateEvaluationData.ObjectiveValueBase;
        end
        optimumCandidateBaseObjectiveInitialValue = optimizationData.IterationData(1).BaseOptimizationData.InitialData.OptimumCandidateEvaluationData.ObjectiveValueBase;
        optimumCandidateBaseObjectiveValue = [optimumCandidateBaseObjectiveInitialValue;optimumCandidateBaseObjectiveValue];
    else
        optimumCandidateEvaluationData = [optimizationData.InitialData.OptimumCandidateEvaluationData;...
            vertcat(optimizationData.IterationData(:).OptimumCandidateEvaluationData)];
        optimumCandidateBaseObjectiveValue = [vertcat(optimumCandidateEvaluationData.ObjectiveValueBase)];
    end
    optimalObjectiveValue = optimumCandidateBaseObjectiveValue(end,:);
    nObjectiveValues = size(optimumCandidateBaseObjectiveValue,2);

    % design variable names
    if(isempty(options.DesignVariableNames))
        for i=1:nDimension
            options.DesignVariableNames{i} = sprintf('$$y_{%d}$$',i);
        end
    end

    % objective value names
    if(isempty(options.ObjectiveValueNames))
        for i=1:nObjectiveValues
            options.ObjectiveValueNames{i} = sprintf('$$y_{%d}$$',i);
        end
    end

    %pairwise plot scatter-plot for multi-objective values
    figure;
    for i = 1:nObjectiveValues
        for j = 1:nObjectiveValues
        subplot(nObjectiveValues, nObjectiveValues, (i-1)*nObjectiveValues + j);
            if i == j
                plot(iterationIndex,optimumCandidateBaseObjectiveValue(:,i), 'LineWidth', options.LineWidth, 'Color', options.Color, ...
                    'Marker', options.Marker, 'LineStyle', options.LineStyle);
                xlabel('Iterations','interpreter','latex','FontSize',options.AxisFontSize)
                ylabel(options.ObjectiveValueNames{i},'interpreter','latex','FontSize',options.AxisFontSize);
                xlim([min(iterationIndex) max(iterationIndex)]);
                if(options.UseDesignSpaceAsAxisLimits)
                    ylim([min(optimumCandidateBaseObjectiveValue(:,i))  max(optimumCandidateBaseObjectiveValue(:,i))]);
                end
                grid(options.ShowGrid);
            else
                scatter(optimumCandidateBaseObjectiveValue(:,j), optimumCandidateBaseObjectiveValue(:,i), options.MarkerSize, optimumCandidateObjectiveValue, 'filled');
                colorbar;
                xlabel(options.ObjectiveValueNames{j},'interpreter','latex','FontSize',options.AxisFontSize);
                ylabel(options.ObjectiveValueNames{i},'interpreter','latex','FontSize',options.AxisFontSize);
                if(options.UseDesignSpaceAsAxisLimits)
                    xlim([min(optimumCandidateBaseObjectiveValue(:,j))  max(optimumCandidateBaseObjectiveValue(:,j))]);
                    ylim([min(optimumCandidateBaseObjectiveValue(:,i))  max(optimumCandidateBaseObjectiveValue(:,i))]);
                end
                grid(options.ShowGrid)
            end
        end
    end
    if(~isempty(options.SaveFolder))
        save_print_figure(gcf,[options.SaveFolder,'PairwiseMatrixObjectiveValues'],options.SaveFigureOptions{:});
        if(options.CloseFigureAfterSaving)
            close(gcf);
        end
    end

    %Plot the objective values of optimal candidate vs iterations
    if(~isempty(options.SaveFolder))
        currentFolder = [options.SaveFolder,'ObjectiveValuesPerIteration/'];
        if(~exist(currentFolder,'dir'))
            mkdir(currentFolder);
        end
    end

    for i = 1:nObjectiveValues
        figure;
        plot(iterationIndex,optimumCandidateBaseObjectiveValue(:,i), 'LineWidth', options.LineWidth, 'Color', options.Color, ...
            'Marker', options.Marker, 'LineStyle', options.LineStyle);
        xlabel('Iterations','interpreter','latex','FontSize',options.AxisFontSize)
        ylabel(options.ObjectiveValueNames{i},'interpreter','latex','FontSize',options.AxisFontSize);
        xlim([min(iterationIndex) max(iterationIndex)]);
        if(options.UseDesignSpaceAsAxisLimits)
            ylim([min(optimumCandidateBaseObjectiveValue(:,i))  max(optimumCandidateBaseObjectiveValue(:,i))]);
        end
        grid(options.ShowGrid)

        if(~isempty(options.SaveFolder))
            filename = sprintf('ObjectiveValue%0*dPerIteration',get_number_digits_integer(nObjectiveValues),i);
            save_print_figure(gcf,[currentFolder,filename],options.SaveFigureOptions{:});
            if(options.CloseFigureAfterSaving)
                close(gcf);
            end
        end
    end

    % normalized objective values per iteration
    optimumCandidateBaseObjectiveValue(optimumCandidateBaseObjectiveValue == inf) = NaN;
    normalizedOptimumCandidateObjectiveValue = (optimumCandidateBaseObjectiveValue-min(optimumCandidateBaseObjectiveValue,[],1))./...
        (max(optimumCandidateBaseObjectiveValue,[],1,'omitnan')-min(optimumCandidateBaseObjectiveValue,[],1,'omitnan'));
    figure;
    plot(iterationIndex,normalizedOptimumCandidateObjectiveValue, 'LineWidth', options.LineWidth, ...
            'Marker', options.Marker, 'LineStyle', options.LineStyle)
    xlabel('Iterations','interpreter','latex','FontSize',options.AxisFontSize)
    ylabel('Normalized Objective Values','interpreter','latex','FontSize',options.AxisFontSize);
    legend(options.ObjectiveValueNames,'interpreter','latex');
    xlim([min(iterationIndex) max(iterationIndex)]);
    ylim([0 1]);
    grid(options.ShowGrid);
    if(~isempty(options.SaveFolder))
        save_print_figure(gcf,[options.SaveFolder,'NormalizedObjectiveValuesPerIteration'],options.SaveFigureOptions{:});
        if(options.CloseFigureAfterSaving)
            close(gcf);
        end
    end

    %error in objective values of optimal candidate vs iterations
    if(~isempty(options.SaveFolder))
        currentFolder = [options.SaveFolder,'ErrorObjectiveValuesPerIteration/'];
        if(~exist(currentFolder,'dir'))
            mkdir(currentFolder);
        end
    end
    for i = 1:nObjectiveValues
        figure;
        plot(iterationIndex, optimumCandidateBaseObjectiveValue(end,i) - optimumCandidateBaseObjectiveValue(:,i), 'LineWidth', options.LineWidth, 'Color', options.Color, ...
            'Marker', options.Marker, 'LineStyle', options.LineStyle);
        xlabel('Iterations','interpreter','latex','FontSize',options.AxisFontSize)
        ylabel(sprintf('Error in %s ($$y^* - y$$)',options.DesignVariableNames{i}),'interpreter','latex','FontSize',options.AxisFontSize);
        xlim([min(iterationIndex) max(iterationIndex)]);
        grid(options.ShowGrid)

        if(~isempty(options.SaveFolder))
            filename = sprintf('ErrorObjectiveValue%0*dPerIteration',get_number_digits_integer(nObjectiveValues),i);
            save_print_figure(gcf,[currentFolder,filename],options.SaveFigureOptions{:});
            if(options.CloseFigureAfterSaving)
                close(gcf);
            end
        end
    end

    % error of normalized objective variable per iteration
    normalizedOptimalObjectiveValue = (optimalObjectiveValue-min(optimumCandidateBaseObjectiveValue,[],1))./...
        (max(optimumCandidateBaseObjectiveValue,[],1)-min(optimumCandidateBaseObjectiveValue,[],1));
    figure;
    plot(iterationIndex,normalizedOptimalObjectiveValue-normalizedOptimumCandidateObjectiveValue, 'LineWidth', options.LineWidth, ...
            'Marker', options.Marker, 'LineStyle', options.LineStyle)
    xlabel('Iterations','interpreter','latex','FontSize',options.AxisFontSize)
    ylabel('Error in Normalized Objective Values','interpreter','latex','FontSize',options.AxisFontSize);
    legend(options.ObjectiveValueNames,'interpreter','latex');
    xlim([min(iterationIndex) max(iterationIndex)]);
    grid(options.ShowGrid);
    if(~isempty(options.SaveFolder))
        save_print_figure(gcf,[options.SaveFolder,'ErrorNormalizedObjectiveValuesPerIteration'],options.SaveFigureOptions{:});
        if(options.CloseFigureAfterSaving)
            close(gcf);
        end
    end
    
end