function plot_point_optimization_results(optimizationData, varargin)
%PLOT_POINT_OPTIMIZATION_RESULTS Visualize optimization trajectory and results for point-based methods.
%
%   PLOT_POINT_OPTIMIZATION_RESULTS(OPTIMIZATIONDATA, ...) generates comprehensive
%   visualizations of optimization progress including design variable evolution,
%   objective function convergence, and 2D/3D contour/surface plots (for n≤3).
%
%   INPUTS:
%       OPTIMIZATIONDATA     - Structure containing optimization problem data and history:
%           .ProblemData     - Optimization problem definition
%           .InitialData     - Initial candidate information
%           .IterationData   - Array of iteration records
%
%   NAME-VALUE PAIR ARGUMENTS:
%       'SaveFolder'         - Path to save figures (default: [])
%       'CloseFigureAfterSaving' - Close figures after saving (default: false)
%       'SaveFigureOptions'  - Cell array of savefig options (default: {})
%       'LineWidth'          - Plot line width (default: 2)
%       'Color'              - Main plot color (default: 'b')
%       'Marker'             - Data point marker (default: 'o')
%       'DesignVariableNames' - Custom names for design variables
%       'ShowContour'        - Show contour plots for 2D problems (default: true)
%       'ShowSurface'        - Show surface plot for 2D problems (default: true)
%       (Complete list of 25+ parameters available in function)
%
%   OUTPUTS:
%       None (generates figures)
%
%   FUNCTIONALITY:
%       The function produces:
%       1. Pairwise variable matrix with iteration history and contour plots
%       2. Individual variable convergence plots
%       3. Normalized design variable evolution
%       4. Objective function convergence
%       5. Error plots for variables and objective
%       6. 2D contour plots with optimization path (for n=2)
%       7. 3D surface plots (for n=2)
%
%   NOTES:
%       - For 2D problems, contour plots use a 100x100 grid by default
%       - All plots use LaTeX interpreters for mathematical notation
%       - Figure saving supports multiple formats via SaveFigureOptions
%       - Axis limits can be locked to design space bounds
%
%   DEPENDENCIES:
%       - save_print_figure.m
%       - get_number_digits_integer.m
%       - evaluate_optimization_objective.m
%       - plot_point_optimization_contour_2d.m
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor, Contributor)
%   Copyright 2025 Ali Abbas Kapadia (Author)
%   Copyright 2025 Gaurav Vaibhav (Author)
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
    addParameter(p, 'LegendFontSize', 12, @(x) isnumeric(x) && isscalar(x));   
    addParameter(p, 'ShowGrid', 'on', @(x) islogical(x) || ischar(x) && any(strcmpi(x, {'on', 'off'})));
    addParameter(p, 'ShowSurface', true);
    addParameter(p, 'ShowSurfaceScatter', false);
    addParameter(p, 'ShowConstraint',true);
    addParameter(p, 'Constraints',[]);
    addParameter(p,'ShowContour', true);
    addParameter(p,'ContourLevels',10);
    addParameter(p,'ContourOptions',{});
    addParameter(p, 'DesignVariableNames',{});
    addParameter(p, 'UseDesignSpaceAsAxisLimits',true);
    addParameter(p, 'PlotContour2DUseDesignSpaceLimit',true);
    addParameter(p,'SubPlotStructure','horizontal',@(x) ischar(x) || isstring(x));
    addParameter(p, 'PlotContour2DOptions',{});
    
    parse(p, varargin{:});
    options = p.Results;

    % set text interpreter
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultAxesFontSize', options.AxisFontSize);
    
    % get problem information / final solution
    nDimension = size(optimizationData.ProblemData(1).DesignSpaceLowerBound,2);
    optimalDesign = optimizationData.IterationData(end).OptimumCandidate;
    optimalObjectiveValue = optimizationData.IterationData(end).OptimumCandidateObjectiveValue;

    % design variable names
    if(isempty(options.DesignVariableNames))
        for i=1:nDimension
            options.DesignVariableNames{i} = sprintf('$$x_{%d}$$',i);
        end
    end

    % extract optimum data
    optimumCandidate = [optimizationData.InitialData.OptimumCandidate;...
        vertcat(optimizationData.IterationData(:).OptimumCandidate)];
    optimumCandidateObjectiveValue = [optimizationData.InitialData.OptimumCandidateObjectiveValue;...
        vertcat(optimizationData.IterationData(:).OptimumCandidateObjectiveValue)];
    iterationIndex = 0:length(optimizationData.IterationData);

    xLimits = [optimizationData.ProblemData(1).DesignSpaceLowerBound(1), optimizationData.ProblemData(1).DesignSpaceUpperBound(1)];
    yLimits = [optimizationData.ProblemData(1).DesignSpaceLowerBound(2), optimizationData.ProblemData(1).DesignSpaceUpperBound(2)];
    if options.ShowContour
        [x,y] = meshgrid(linspace(xLimits(1),xLimits(2),100),linspace(yLimits(1),yLimits(2),100));
        fullGrid = [x(:),y(:)];
        z = evaluate_optimization_objective(optimizationData.ProblemData(1).ObjectiveFunction,fullGrid);
        z = reshape(z, size(x));
        % Set colormap and global color axis
        colormap('jet');
        cmin = min(z(:));
        cmax = max(z(:));
        clim([cmin cmax]);                                                 % Same color range for all subplots
    end
  
    % pairwise matrix scatter-plot for all dimensions
    ax = gobjects(nDimension, nDimension);
    x_labels = cell(nDimension, nDimension);
    y_labels = cell(nDimension, nDimension);
    for i = 1:nDimension
        for j = 1:nDimension
            ax(i,j) = subplot(nDimension, nDimension, (i-1)*nDimension + j);
            if i == j
                plot(iterationIndex,optimumCandidate(:,i), 'LineWidth', options.LineWidth, 'Color', options.Color, ...
                    'Marker', options.Marker, 'LineStyle', options.LineStyle);
                xlabel('Iteration Index (k)','interpreter','latex','FontSize',options.AxisFontSize);
                ylabel(options.DesignVariableNames{i},'interpreter','latex','FontSize',options.AxisFontSize);
                xlim([min(iterationIndex) max(iterationIndex)]);
                xticks(min(iterationIndex):max(iterationIndex));
                xticklabels(string(min(iterationIndex):max(iterationIndex)));  % Force integer labelsbels
                if(options.UseDesignSpaceAsAxisLimits)
                    ylim([optimizationData.ProblemData(1).DesignSpaceLowerBound(i) optimizationData.ProblemData(1).DesignSpaceUpperBound(i)]);
                end
                grid(options.ShowGrid);
            else
                hold on;
                scatter(optimumCandidate(:,j), optimumCandidate(:,i), options.MarkerSize, optimumCandidateObjectiveValue, 'filled');
                colorbar;
                if options.ShowContour
                    if(i==1)
                        contour(y,x,z,options.ContourLevels,options.ContourOptions{:});
                        colormap('jet');
                    else
                        contour(x,y,z,options.ContourLevels,options.ContourOptions{:});
                        colormap('jet');
                    end
                end
                xlabel(options.DesignVariableNames{j},'interpreter','latex','FontSize',options.AxisFontSize);
                ylabel(options.DesignVariableNames{i},'interpreter','latex','FontSize',options.AxisFontSize);
                if(options.UseDesignSpaceAsAxisLimits)
                    xlim([optimizationData.ProblemData(1).DesignSpaceLowerBound(j) optimizationData.ProblemData(1).DesignSpaceUpperBound(j)]);
                    ylim([optimizationData.ProblemData(1).DesignSpaceLowerBound(i) optimizationData.ProblemData(1).DesignSpaceUpperBound(i)]);
                end
                grid(options.ShowGrid)
            end
        end
    end
    % ----------------------------
    % Link axes
    % ----------------------------
    % first, link x-axes of all diagonals (Iterations)
    for i = 1:nDimension-1
        linkaxes([ax(i,i), ax(i+1,i+1)], 'x');
    end
    % now, link y-axis of everything in the same row (same variable)
    for i=1:nDimension
        for j=i+1:nDimension
            linkaxes([ax(i,j-1), ax(i,j)], 'y');
        end
    end
    % now do cross-linking (link y-axis of the diagonal plot in a column to x-axis of all other entries in the same column)
    for j=1:nDimension
        for i=1:nDimension
            if(i~=j)
                addlistener(ax(j,j), 'YLim', 'PostSet', @(src, evnt) set(ax(i,j), 'XLim', get(ax(j,j), 'YLim')));
            end
        end
    end
    if(~isempty(options.SaveFolder))
        save_print_figure(gcf,[options.SaveFolder,'PairwiseMatrix'],options.SaveFigureOptions{:});
        if(options.CloseFigureAfterSaving)
            close(gcf);
        end
    end

    
    %Plot the design variables values of optimal candidate vs iterations
    if(~isempty(options.SaveFolder))
        currentFolder = [options.SaveFolder,'DesignVariablePerIteration/'];
        if(~exist(currentFolder,'dir'))
            mkdir(currentFolder);
        end
    end
    for i = 1:nDimension
        figure;
        set(gcf, 'Color', 'w');  % Set the figure background to white
        plot(iterationIndex,optimumCandidate(:,i), 'LineWidth', options.LineWidth, 'Color', options.Color, ...
            'Marker', options.Marker, 'LineStyle', options.LineStyle);
        xlabel('Iteration Index (k)','interpreter','latex','FontSize',options.AxisFontSize)
        ylabel(options.DesignVariableNames{i},'interpreter','latex','FontSize',options.AxisFontSize);
        xlim([min(iterationIndex) max(iterationIndex)]);
        xticks(min(iterationIndex):max(iterationIndex));
        xticklabels(string(min(iterationIndex):max(iterationIndex)));  % Force integer labels

        if(options.UseDesignSpaceAsAxisLimits)
            ylim([optimizationData.ProblemData(1).DesignSpaceLowerBound(i) optimizationData.ProblemData(1).DesignSpaceUpperBound(i)]);
        end
        grid(options.ShowGrid)

        if(~isempty(options.SaveFolder))
            filename = sprintf('DesignVariable%0*dPerIteration',get_number_digits_integer(nDimension),i);
            save_print_figure(gcf,[currentFolder,filename],options.SaveFigureOptions{:});
            if(options.CloseFigureAfterSaving)
                close(gcf);
            end
        end
    end

    % normalized design variables per iteration
    normalizedOptimumCandidate = (optimumCandidate-optimizationData.ProblemData(1).DesignSpaceLowerBound)./...
        (optimizationData.ProblemData(1).DesignSpaceUpperBound-optimizationData.ProblemData(1).DesignSpaceLowerBound);
    figure;
    set(gcf, 'Color', 'w');  % Set the figure background to white
    if i == 1
        plot(iterationIndex,normalizedOptimumCandidate, 'LineWidth', options.LineWidth,'Color', options.Color, ...
                'Marker', options.Marker, 'LineStyle', options.LineStyle);
    else
        plot(iterationIndex,normalizedOptimumCandidate, 'LineWidth', options.LineWidth, ...
        'Marker', options.Marker, 'LineStyle', options.LineStyle);
    end
    xlabel('Iteration Index (k)','interpreter','latex','FontSize',options.AxisFontSize);
    ylabel('Normalized Design Variables','interpreter','latex','FontSize',options.AxisFontSize);
    legend(options.DesignVariableNames,'interpreter','latex','FontSize',options.LegendFontSize);
    xlim([min(iterationIndex) max(iterationIndex)]);
    xticks(min(iterationIndex):max(iterationIndex));
    xticklabels(string(min(iterationIndex):max(iterationIndex)));  % Force integer labels

    ylim([0 1]);
    grid(options.ShowGrid);
    if(~isempty(options.SaveFolder))
        save_print_figure(gcf,[options.SaveFolder,'NormalizedDesignVariablesPerIteration'],options.SaveFigureOptions{:});
        if(options.CloseFigureAfterSaving)
            close(gcf);
        end
    end

    %Plot the error in design variables values of optimal candidate vs iterations
    if(~isempty(options.SaveFolder))
        currentFolder = [options.SaveFolder,'ErrorDesignVariablePerIteration/'];
        if(~exist(currentFolder,'dir'))
            mkdir(currentFolder);
        end
    end
    for i = 1:nDimension
        figure;
        set(gcf, 'Color', 'w');  % Set the figure background to white
        plot(iterationIndex, optimalDesign(i) - optimumCandidate(:,i), 'LineWidth', options.LineWidth, 'Color', options.Color, ...
            'Marker', options.Marker, 'LineStyle', options.LineStyle);
        xlabel('Iteration Index (k)','interpreter','latex','FontSize',options.AxisFontSize)
        ylabel(sprintf('Error in %s ($$x^* - x$$)',options.DesignVariableNames{i}),'interpreter','latex','FontSize',options.AxisFontSize);
        xlim([min(iterationIndex) max(iterationIndex)]);
        xticks(min(iterationIndex):max(iterationIndex));
        xticklabels(string(min(iterationIndex):max(iterationIndex)));  % Force integer labels
        grid(options.ShowGrid)

        if(~isempty(options.SaveFolder))
            filename = sprintf('ErrorDesignVariable%0*dPerIteration',get_number_digits_integer(nDimension),i);
            save_print_figure(gcf,[currentFolder,filename],options.SaveFigureOptions{:});
            if(options.CloseFigureAfterSaving)
                close(gcf);
            end
        end
    end

    % error normalized design variable per iteration
    normalizedOptimalDesign = (optimalDesign-optimizationData.ProblemData(1).DesignSpaceLowerBound)./...
        (optimizationData.ProblemData(1).DesignSpaceUpperBound-optimizationData.ProblemData(1).DesignSpaceLowerBound);
    figure;
    set(gcf, 'Color', 'w');  % Set the figure background to white
    if i == 1
        plot(iterationIndex,normalizedOptimumCandidate, 'LineWidth', options.LineWidth,'Color', options.Color, ...
                'Marker', options.Marker, 'LineStyle', options.LineStyle);
    else
        plot(iterationIndex,normalizedOptimumCandidate, 'LineWidth', options.LineWidth, ...
        'Marker', options.Marker, 'LineStyle', options.LineStyle);
    end
    xlabel('Iteration Index (k)','interpreter','latex','FontSize',options.AxisFontSize);
    ylabel('Error in Normalized Design Variables','interpreter','latex','FontSize',options.AxisFontSize);
    legend(options.DesignVariableNames,'interpreter','latex','FontSize',options.LegendFontSize);
    xlim([min(iterationIndex) max(iterationIndex)]);
    xticks(min(iterationIndex):max(iterationIndex));
    xticklabels(string(min(iterationIndex):max(iterationIndex)));  % Force integer labels
    grid(options.ShowGrid);
    if(~isempty(options.SaveFolder))
        save_print_figure(gcf,[options.SaveFolder,'ErrorNormalizedDesignVariablesPerIteration'],options.SaveFigureOptions{:});
        if(options.CloseFigureAfterSaving)
            close(gcf);
        end
    end
    
    %Plot the Quantity of interest value of optimal candidate vs iterations
    figure;
    set(gcf, 'Color', 'w');  % Set the figure background to white
    plot(iterationIndex,optimumCandidateObjectiveValue, 'LineWidth', options.LineWidth, 'Color', options.Color, ...
        'Marker', options.Marker, 'LineStyle', options.LineStyle);
    xlabel('Iteration Index (k)','interpreter','latex','FontSize',options.AxisFontSize);
    ylabel('Objective Value ($f(\mathbf{x})$)','interpreter','latex','FontSize',options.AxisFontSize);
    xlim([min(iterationIndex) max(iterationIndex)]);
    xticks(min(iterationIndex):max(iterationIndex));
    xticklabels(string(min(iterationIndex):max(iterationIndex)));  % Force integer labels
    grid(options.ShowGrid);
    if(~isempty(options.SaveFolder))
        save_print_figure(gcf,[options.SaveFolder,'ObjectiveValuePerIteration'],options.SaveFigureOptions{:});
        if(options.CloseFigureAfterSaving)
            close(gcf);
        end
    end
    
    %Plot the error in Quantity of interest value of optimal candidate vs iterations
    figure;
    set(gcf, 'Color', 'w');  % Set the figure background to white
    plot(iterationIndex,optimalObjectiveValue - optimumCandidateObjectiveValue, 'LineWidth', options.LineWidth, 'Color', options.Color, ...
            'Marker', options.Marker, 'LineStyle', options.LineStyle);
    xlabel('Iteration Index (k)','interpreter','latex','FontSize',options.AxisFontSize);
    ylabel('Error in Objective Function ($$f^* - f$$)','interpreter','latex','FontSize',options.AxisFontSize);
    xlim([min(iterationIndex) max(iterationIndex)]);
    xticks(min(iterationIndex):max(iterationIndex));
    xticklabels(string(min(iterationIndex):max(iterationIndex)));  % Force integer labels
    grid(options.ShowGrid);
    if(~isempty(options.SaveFolder))
        save_print_figure(gcf,[options.SaveFolder,'ErrorObjectiveValuePerIteration'],options.SaveFigureOptions{:});
        if(options.CloseFigureAfterSaving)
            close(gcf);
        end
    end
    
    % defining rows, cols, and subPlotIndex function
    switch lower(options.SubPlotStructure)
        case 'vertical'
            nRows = nDimension; % Stack all plots vertically
            nCols = 3;
            subPlotIndex = @(i, j) (i - 1) * 3 + j;
        case 'horizontal'
            nRows = 3;
            nCols = nDimension; % Arrange all plots horizontally
            subPlotIndex = @(i, j) (j - 1) * nDimension + i;
        case 'separate'
            nRows = 3;
            nCols = 1;
            subPlotIndex = @(i,j) j;
    end
    
    % Create a tiled layout with consistent spacing
    tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');  % or 'none' for tighter layout
    
    ax = gobjects(nDimension, 3); % Preallocate
    ylimData = zeros(nRows, 2, 2); % To store ylim of each subplot (only 2 columns per row)
    
    % Loop for plotting
    for i = 1:nDimension
        % Plot 1: Design Variable
        ax(i,1) = nexttile(subPlotIndex(i,1));
        plot(iterationIndex, optimumCandidate(:,i), 'LineWidth', options.LineWidth, 'Color', options.Color, ...
            'Marker', options.Marker, 'LineStyle', options.LineStyle);
        set(gcf, 'Color', 'w');  % Set the figure background to white
        xlabel('Iteration Index (k)', 'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
        ylabel(options.DesignVariableNames{i}, 'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
        xlim([min(iterationIndex) max(iterationIndex)]);
        xticks(min(iterationIndex):max(iterationIndex));
        xticklabels(string(min(iterationIndex):max(iterationIndex)));
        if options.UseDesignSpaceAsAxisLimits
            ylim([optimizationData.ProblemData(1).DesignSpaceLowerBound(i), ...
                  optimizationData.ProblemData(1).DesignSpaceUpperBound(i)]);
        end
        ylimData(i,1,:) = ylim; % Save y-limits for the first plot
        grid(options.ShowGrid);
    
        % Plot 2: Error
        ax(i,2) = nexttile(subPlotIndex(i,2));
        plot(iterationIndex, optimalDesign(i) - optimumCandidate(:,i), 'LineWidth', options.LineWidth, 'Color', options.Color, ...
            'Marker', options.Marker, 'LineStyle', options.LineStyle);
        set(gcf, 'Color', 'w');  % Set the figure background to white
        xlabel('Iteration Index (k)', 'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
        ylabel(sprintf('Error in %s ($$x^* - x$$)', options.DesignVariableNames{i}), ...
            'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
        xlim([min(iterationIndex) max(iterationIndex)]);
        xticks(min(iterationIndex):max(iterationIndex));
        xticklabels(string(min(iterationIndex):max(iterationIndex)));
        ylimData(i,2,:) = ylim; % Save y-limits for the second plot
        grid(options.ShowGrid);
    
        % Plot 3: Normalized
        ax(i,3) = nexttile(subPlotIndex(i,3));
        plot(iterationIndex, normalizedOptimumCandidate(:,i), 'LineWidth', options.LineWidth, 'Color', options.Color, ...
            'Marker', options.Marker, 'LineStyle', options.LineStyle);
        set(gcf, 'Color', 'w');  % Set the figure background to white
        xlabel('Iteration Index (k)', 'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
        ylabel(sprintf('Normalized %s', options.DesignVariableNames{i}), 'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
        xlim([min(iterationIndex) max(iterationIndex)]);
        xticks(min(iterationIndex):max(iterationIndex));
        xticklabels(string(min(iterationIndex):max(iterationIndex)));
        ylim([0 1]);  % fixed for normalized plots
        ylimData(i,3,:) = [0 1]; % Store fixed limits for normalized plot
        grid(options.ShowGrid);
    end
    
    % Adjust Y-axis to match grid divisions within each row (same y-ticks for both plots in a row)
    for row = 1:nDimension
        % Get min and max y-limits for the row (across the two first plots)
        yMin = min([ylimData(row, 1, 1), ylimData(row, 2, 1)]);
        yMax = max([ylimData(row, 1, 2), ylimData(row, 2, 2)]);
        
        % Get the number of ticks from the first plot in the row (to maintain consistent grid divisions)
        numTicks = numel(get(ax(row, 1), 'YTick'));
        
        % Create shared ticks for the row
        sharedTicks = linspace(yMin, yMax, numTicks);
        
        % Apply the same Y-limits and ticks to both plots in the row
        for col = 1:2
            set(ax(row, col), 'YLim', [yMin yMax], 'YTick', sharedTicks);
        end
    end
    
    % Synchronize x-axes across all plots
    linkaxes(ax(:), 'x');

    for k = 1:numel(ax)
        yl = ax(k).YLabel;
        pos = yl.Position;
        pos(1) = -0.25;  % Adjust horizontal position
        yl.Position = pos;
    end


    %ContourPlot for dimension=2
    if nDimension == 2
        figureHandle = figure;

        if(options.PlotContour2DUseDesignSpaceLimit)
            contourDesignSpace = [optimizationData.ProblemData(1).DesignSpaceLowerBound;optimizationData.ProblemData(1).DesignSpaceUpperBound];
        else
            contourDesignSpace = [min(optimumCandidate,[],1);max(optimumCandidate,[],1)];
        end
        plot_point_optimization_contour_2d(figureHandle, optimizationData.ProblemData(1).ObjectiveFunction, contourDesignSpace, ...
            'DesignVariableNames',options.DesignVariableNames,'AxisFontSize',options.AxisFontSize,...
            options.PlotContour2DOptions{:});

        % Gradient -> normalized
        gradientDirection = optimizationData.ProblemData(1).DesignSpaceLowerBound - optimizationData.IterationData(1).OptimumCandidate;  % Ensure column
        gradientDescent = -1 * (gradientDirection / norm(gradientDirection + eps));  % Normalize with epsilon for safety

        gradientStartPoint = optimizationData.ProblemData(1).DesignSpaceLowerBound + 0.1;
        gradientEndPoint = gradientStartPoint + 0.25*gradientDescent;
        
        % Plot arrow
        hold on;
        h_quiver = quiver(gradientStartPoint(1), gradientStartPoint(2), ...
                       gradientEndPoint(1) - gradientStartPoint(1), ...
                       gradientEndPoint(2) - gradientStartPoint(2), ...
                       0, 'k', 'LineWidth', 2, 'MaxHeadSize', 2, ...
                       'DisplayName', 'Gradient');

        set(gcf, 'Color', 'w');  % Set the figure background to white
        %scatter(optimumCandidate(:,j), optimumCandidate(:,i), options.MarkerSize, optimumCandidateObjectiveValue, 'filled');
        if(options.UseDesignSpaceAsAxisLimits)
            xlim([optimizationData.ProblemData(1).DesignSpaceLowerBound(1) optimizationData.ProblemData(1).DesignSpaceUpperBound(1)]);
            ylim([optimizationData.ProblemData(1).DesignSpaceLowerBound(2) optimizationData.ProblemData(1).DesignSpaceUpperBound(2)]);
        end
        grid(options.ShowGrid)
        legend(h_quiver, 'Location', 'northwest');

        if(~isempty(options.SaveFolder))
            save_print_figure(gcf,[options.SaveFolder,'ContourPlot2D'],options.SaveFigureOptions{:});
            if(options.CloseFigureAfterSaving)
                close(gcf);
            end
        end

        if options.ShowSurface
            figureHandle = figure;
            ax = axes(figureHandle);
            surf(ax,x,y,z,'EdgeColor','none');
            hold on;
            colormap(ax,'jet');
            colorbar(ax);
            xlabel(options.DesignVariableNames{1},'interpreter','latex','FontSize',options.AxisFontSize);
            ylabel(options.DesignVariableNames{2},'interpreter','latex','FontSize',options.AxisFontSize);
            zlabel('Objective Value ($f(\mathbf{x})$)','interpreter','latex','FontSize',options.AxisFontSize);
            
            if options.ShowSurfaceScatter
                scatter3(optimumCandidate(:,1), optimumCandidate(:,2),optimumCandidateObjectiveValue,options.MarkerSize,'w','filled');
            end

            if(options.UseDesignSpaceAsAxisLimits)
                xlim([optimizationData.ProblemData(1).DesignSpaceLowerBound(1) optimizationData.ProblemData(1).DesignSpaceUpperBound(1)]);
                ylim([optimizationData.ProblemData(1).DesignSpaceLowerBound(2) optimizationData.ProblemData(1).DesignSpaceUpperBound(2)]);
            end
            grid(options.ShowGrid);
            view(3);
            if(~isempty(options.SaveFolder))
                save_print_figure(gcf,[options.SaveFolder,'SurfacePlot3D'],options.SaveFigureOptions{:});
                if(options.CloseFigureAfterSaving)
                    close(gcf);
                end
            end
        end

        
    end
end

