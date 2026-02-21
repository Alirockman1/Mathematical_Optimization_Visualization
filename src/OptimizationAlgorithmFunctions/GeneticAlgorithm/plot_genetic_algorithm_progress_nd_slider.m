function plot_genetic_algorithm_progress_nd_slider(optimizationData, varargin)
%PLOT_GENETIC_ALGORITHM_PROGRESS_ND_SLIDER Visualize GA progress in nD design space using 2D projections.
%
%   PLOT_GENETIC_ALGORITHM_PROGRESS_ND_SLIDER(DATA, VARIABLE_NAMES, ...
%       OBJECTIVE_NAME, OPTIONS) creates interactive figures for all 
%       variable pairs in an n-dimensional genetic algorithm, showing 
%       the evolution of different GA steps (Current, Selection, 
%       Crossover, Mutation) across iterations using subplots, contour 
%       overlays, and a slider for iteration control.
%
%   INPUTS:
%       DATA                - Struct with GA population data:
%           .current        - [N x D x I] array of current populations
%           .selection      - [N x D x I] array of selected parents
%           .crossover      - [N x D x I] array of crossover offspring
%           .mutation       - [N x D x I] array of mutated candidates
%           .objective      - [N x 1 x I] array of objective values (optional)
%
%       VARIABLE_NAMES      - Cell array of strings with names of each design variable
%       OBJECTIVE_NAME      - String for labeling the objective (e.g., 'f(x)')
%
%   NAME-VALUE PAIR ARGUMENTS (inside OPTIONS struct):
%       'ContourFunction'   - Function handle @(x, y) for background contours (optional)
%       'ContourResolution' - Integer resolution for contour grid (default: 100)
%       'MarkerSize'        - Marker size for scatter plots (default: 25)
%       'LineWidth'         - Border width for markers (default: 1.5)
%       'ShowLegend'        - Logical flag to show legend (default: true)
%       'Transparency'      - Alpha transparency for markers (default: 0.7)
%       'ColorMap'          - Colormap name or [N x 3] RGB array (default: 'lines')
%
%   OUTPUTS:
%       None. Interactive figures are displayed with one figure per 
%       variable pair, containing four subplots and a UI slider to 
%       navigate iterations.
%
%   FUNCTIONALITY:
%       1. Generates 2D scatter plots for each GA step and variable pair
%       2. Supports visual overlays of user-defined objective contours
%       3. Efficiently updates plots with iteration slider using handle graphics
%       4. Fully customizable visual style via 'options' structure
%       5. Compatible with high-dimensional GA visualization workflows
%
%   NOTES:
%       - Requires MATLAB R2018b or later for uislider and tiledlayout support
%       - Each figure represents a variable pair (i,j), resulting in D*(D-1)/2 plots
%       - Designed for offline inspection or visual reporting of GA behavior
%       - Optimized for visual clarity with minimal redraw overhead
%
%   SEE ALSO:
%       plot_objective_value_progress_nd, plot_differential_evolution_progress_2d_slider
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor, Contributor)
%   Copyright 2025 Ali Abbas Kapadia (Author) 
%   Copyright 2025 Gaurav Vaibhav (Author)
%   SPDX-License-Identifier: Apache-2.0

    % Parse input parameters
    parser = inputParser;
    parser.addParameter('SaveFolder', []);
    parser.addParameter('SaveFigureOptions', {});
    parser.addParameter('CreateVideoGifIterations', true);
    parser.addParameter('FramesPerSecond', 5);
    parser.addParameter('ImwriteOptions', {});
    parser.addParameter('VideoWriterProfile', 'Motion JPEG 2000');
    parser.addParameter('VideoWriterOptions', {'LosslessCompression', true});
    parser.addParameter('ShowContour', true);
    parser.addParameter('ContourLevels', 10);
    parser.addParameter('ContourOptions', {});
    parser.addParameter('IncludeLegend', true);
    parser.addParameter('IncludeArrow', true);
    parser.addParameter('LegendOptions', {});
    parser.addParameter('ShowGrid', 'on', @(x) islogical(x) || ischar(x) && any(strcmpi(x, {'on', 'off'})));
    parser.addParameter('MarkerEdgeColor', 'b');
    parser.addParameter('Marker', 'o');
    parser.addParameter('LineStyle', '-');
    parser.addParameter('LineType', 'line');
    parser.addParameter('MarkerFaceColor', 'auto');
    parser.addParameter('MarkerSize', 36);
    parser.addParameter('AxisFontSize', 12);
    parser.addParameter('LegendFontSize', 12);
    parser.addParameter('DesignVariableNames', {});
    parser.addParameter('UseDesignSpaceAsAxisLimits', true);
    parser.addParameter('PlotContour2DUseDesignSpaceLimit', true);
    parser.addParameter('PlotContour2DOptions', {});
    parser.parse(varargin{:});
    options = parser.Results;

    % Set up design variable names
    if isempty(options.DesignVariableNames)
        for i = 1:length(optimizationData.InitialData.OptimumCandidate)
            options.DesignVariableNames{i} = sprintf('$$x_{%d}$$', i);
        end
    end

    defaultVideoWriterOptions = {'FrameRate',options.FramesPerSecond};
    [~,videoWriterOptions] = merge_name_value_pair_argument(defaultVideoWriterOptions,options.VideoWriterOptions);
    defaultImwriteOptions = {};
    [~,imwriteOptions] = merge_name_value_pair_argument(defaultImwriteOptions,options.ImwriteOptions);
    frameDelay = 1/options.FramesPerSecond;
    defaultLegendOptions = {'location','southeast','FontSize', 6};
    [~,legendOptions] = merge_name_value_pair_argument(defaultLegendOptions,options.LegendOptions);

    % Define the color sequences for each GA step (1 to 4)
    colorSequences = {
        {'b'},             % Step 1: Current Population
        {'b', 'r'},        % Step 2: Selection (prev + selected)
        {'r', 'k'},        % Step 3: Crossover (prev + crossover)
        {'k', 'g'}         % Step 4: Mutation (prev + mutated)
    };

    % Problem info
    designVariables = size(optimizationData.InitialData.PopulationDesignInitial, 2);
    varPairs = nchoosek(1:designVariables, 2);
    nPairs = size(varPairs, 1);
    nIter = length(optimizationData.IterationData);
    subplotTitles = {'Current Population', 'Selection Step', 'CrossOver Step', 'Mutation Step'};

% --------------------------------------------------------------
% Define Plot layout
% --------------------------------------------------------------
    % Create figure
    figureHandle = figure('Name', 'Genetic algorithm progress', 'NumberTitle', 'off');
    set(gcf, 'Color', 'w');    

    % Preallocate axes handles: rows = pairs, cols = 4 GA steps
    axesHandles = gobjects(nPairs, 4);
  
    if designVariables == 2
        % Create subplots for the 4 steps
        ax1 = subplot(2, 2, 1); hold on;
        ax2 = subplot(2, 2, 2); hold on;
        ax3 = subplot(2, 2, 3); hold on;
        ax4 = subplot(2, 2, 4); hold on;

        axesHandles = [ax1,ax2,ax3,ax4];

        % Set axis limits
        xLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(1), optimizationData.ProblemData.DesignSpaceUpperBound(1)];
        yLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(2), optimizationData.ProblemData.DesignSpaceUpperBound(2)];

        % Titles, labels, limits for each subplot
        for i = 1:4
            subplot(2, 2, i);
            title(subplotTitles{i});
            xlabel('x1','interpreter','latex','FontSize',options.AxisFontSize);
            ylabel('x2','interpreter','latex','FontSize',options.AxisFontSize);
            xlim(xLimits);
            ylim(yLimits);
            grid(options.ShowGrid);
        end

        % Plot 2D contour slices
        if options.ShowContour
            % 1. Evaluate objective function grid
            [x, y] = meshgrid(linspace(xLimits(1), xLimits(2), 100), ...
                              linspace(yLimits(1), yLimits(2), 100));
            fullGrid = [x(:), y(:)];
            z = evaluate_optimization_objective(optimizationData.ProblemData.ObjectiveFunction, fullGrid);
            z = reshape(z, size(x));
        
            % 2. Create invisible figure + axis to generate contour objects
            tempFig = figure('Visible', 'off');
            tempAx = axes(tempFig);
            [~, hContourLines] = contour(tempAx, x, y, z, options.ContourLevels, options.ContourOptions{:});
            colormap('jet'); % Apply colormap globally if needed
            cbar = colorbar;
            cbar.Label.String = 'Objective Value';
        
            % 3. For each plot in a row, copy contour lines
            for axesRowIdx = 1:4
                ax = axesHandles(1,axesRowIdx);
                copyobj(hContourLines, ax);
            end
        
            % 5. Clean up
            close(tempFig);  % Close temporary figure
        end

    elseif designVariables > 2
        % Create subplots in grid: nPairs rows × 4 cols (Current, Selection, Crossover, Mutation)
        for pairIdx = 1:nPairs
            dimension1 = varPairs(pairIdx, 1);
            dimension2 = varPairs(pairIdx, 2);
            
            % Set axis limits
            xLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(dimension1), optimizationData.ProblemData.DesignSpaceUpperBound(dimension1)];
            yLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(dimension2), optimizationData.ProblemData.DesignSpaceUpperBound(dimension2)];

            for rowIdx = 1:4
                ax = subplot(nPairs, 4, (pairIdx - 1) * 4 + rowIdx); hold(ax, 'on');
                title(subplotTitles{rowIdx});
                axesHandles(pairIdx, rowIdx) = ax;
                xlabel(sprintf('x%d', varPairs(pairIdx,1)),'interpreter','latex','FontSize',options.AxisFontSize);
                ylabel(sprintf('x%d', varPairs(pairIdx,2)),'interpreter','latex','FontSize',options.AxisFontSize);
                xlim(xLimits);
                ylim(yLimits);
                grid(options.ShowGrid);
            end

            % Plot 2D contour slices
            if options.ShowContour
                % 1. Evaluate objective function grid
                [x, y] = meshgrid(linspace(xLimits(1), xLimits(2), 100), ...
                                  linspace(yLimits(1), yLimits(2), 100));
                fullGrid = [x(:), y(:)];
                z = evaluate_optimization_objective(optimizationData.ProblemData.ObjectiveFunction, fullGrid);
                z = reshape(z, size(x));
            
                % 2. Create invisible figure + axis to generate contour objects
                tempFig = figure('Visible', 'off');
                tempAx = axes(tempFig);
                [~, hContourLines] = contour(tempAx, x, y, z, options.ContourLevels, options.ContourOptions{:});
                colormap('jet'); % Apply colormap globally if needed
                cbar = colorbar;
                cbar.Label.String = 'Objective Value';
            
                % 3. For each plot in a row, copy contour lines
                for axesRowIdx = 1:4
                    ax = axesHandles(pairIdx,axesRowIdx);
                    copyobj(hContourLines, ax);
                end
            
                % 5. Clean up
                close(tempFig);  % Close temporary figure
            end
        end
    else
        error('At least two design variables are required for 2D plotting.');
    end

    % Create empty plot handles to update later
    plotHandles = cell(nPairs, 4);  % Rows: variable pairs, Cols: GA steps

    for pairIdx = 1:nPairs
        for step = 1:4
            ax = axesHandles(pairIdx, step);
            nPlots = numel(colorSequences{step});
            plotHandles{pairIdx, step} = gobjects(1, nPlots);
    
            for i = 1:nPlots
                plotHandles{pairIdx, step}(i) = plot(ax, NaN, NaN, '.', ...
                    'MarkerSize', 10, ...
                    'Color', colorSequences{step}{i});
            end
        end
    end

% --------------------------------------------------------------
% Define Slider UI element
% --------------------------------------------------------------
    slider = uicontrol('Style', 'slider', ...
        'Min', 1, 'Max', nIter, ...
        'Value', 1, 'SliderStep', [1/(nIter-1), 1/(nIter-1)], ...
        'Position', [100, 10, 400, 20], ...
        'Callback', @(src, event) updatePlot(round(src.Value)));

% --------------------------------------------------------------
% Define Graph Update Function
% --------------------------------------------------------------
    function updatePlot(iter)
        if iter == 1
            populationPrevious = optimizationData.InitialData.PopulationDesignInitial;
        else
            populationPrevious = optimizationData.IterationData(iter-1).MutationPopulation;
        end
        selectionIndex = optimizationData.IterationData(iter).SelectionIndex;
        isMutate = optimizationData.IterationData(iter).MutationIsMutated;
        selectedPopulation = optimizationData.IterationData(iter).SelectionPopulation;
        crossOverPopulation = optimizationData.IterationData(iter).CrossoverPopulation;
        mutatedPopulation = optimizationData.IterationData(iter).MutationPopulation;
    
        sgtitle(sprintf('Genetic Algorithm Iteration %d, Population Size %d', iter, optimizationData.ProblemData.Options.PopulationSize)); % Overarching title
    
        for plotIter = 1:nPairs
            i1 = varPairs(plotIter, 1);
            i2 = varPairs(plotIter, 2);
    
            % Current Population subplot
            ax1 = axesHandles(plotIter, 1);
            set(plotHandles{plotIter,1}(1), 'XData', populationPrevious(:, i1), 'YData', populationPrevious(:, i2));
            if options.IncludeLegend
                legend1 = legend(ax1, plotHandles{plotIter,1}(1), {'Current Population'}, legendOptions{:});
                set(legend1, 'Box', 'on', 'Color', [0.9 0.9 0.9], 'EdgeColor', 'k');
            end
    
            % Selection Step subplot
            ax2 = axesHandles(plotIter, 2);
            set(plotHandles{plotIter,2}(1), 'XData', populationPrevious(:, i1), 'YData', populationPrevious(:, i2));
            set(plotHandles{plotIter,2}(2), 'XData', selectedPopulation(selectionIndex, i1), 'YData', selectedPopulation(selectionIndex, i2));
            if options.IncludeLegend
                legend2 = legend(ax2, [plotHandles{plotIter,2}(1), plotHandles{plotIter,2}(2)], ...
                    {'Current Population', 'Selected Population'}, legendOptions{:});
                set(legend2, 'Box', 'on', 'Color', [0.9 0.9 0.9], 'EdgeColor', 'k');
            end
    
            % Crossover Step subplot
            ax3 = axesHandles(plotIter, 3);
            set(plotHandles{plotIter,3}(1), 'XData', selectedPopulation(:, i1), 'YData', selectedPopulation(:, i2));
            set(plotHandles{plotIter,3}(2), 'XData', crossOverPopulation(:, i1), 'YData', crossOverPopulation(:, i2));
            if options.IncludeLegend
                legend3 = legend(ax3, [plotHandles{plotIter,3}(1), plotHandles{plotIter,3}(2)], ...
                    {'Selected Population', 'Crossover Population'}, legendOptions{:});
                set(legend3, 'Box', 'on', 'Color', [0.9 0.9 0.9], 'EdgeColor', 'k');
            end
    
            % Mutation Step subplot
            ax4 = axesHandles(plotIter, 4);
            set(plotHandles{plotIter,4}(1), 'XData', crossOverPopulation(~isMutate, i1), 'YData', crossOverPopulation(~isMutate, i2));
            set(plotHandles{plotIter,4}(2), 'XData', mutatedPopulation(isMutate, i1), 'YData', mutatedPopulation(isMutate, i2));
            if options.IncludeLegend
                legend4 = legend(ax4, [plotHandles{plotIter,4}(1), plotHandles{plotIter,4}(2)], ...
                    {'Chosen from Crossover', 'Mutated Population'}, legendOptions{:});
                set(legend4, 'Box', 'on', 'Color', [0.9 0.9 0.9], 'EdgeColor', 'k');
            end
        end
    
        drawnow;
    end

% --------------------------------------------------------------
% Plot Graph
% --------------------------------------------------------------
    updatePlot(1);

% --------------------------------------------------------------
% Record Video
% --------------------------------------------------------------
    % Video Creation
    if options.CreateVideoGifIterations
        filename = [options.SaveFolder, 'IterationProgression'];
        videoHandle = VideoWriter(filename, options.VideoWriterProfile);
        for i = 1:2:length(options.VideoWriterOptions)
            videoHandle.(options.VideoWriterOptions{i}) = options.VideoWriterOptions{i+1};
        end

        open(videoHandle);
        
        % Iteration loop
        for playback = 1:nIter
            slider.Value = playback;
            updatePlot(playback);
            currentFrame = getframe(figureHandle);
            writeVideo(videoHandle, currentFrame);

            % GIF Creation
            % Convert to indexed image for GIF
            [A, map] = rgb2ind(frame2im(currentFrame), 256);
            gifPath = [filename, '.gif'];
    
            if playback == 1
                % First frame: create new GIF
                imwrite(A, map, gifPath, 'gif', 'LoopCount', inf, ...
                    'DelayTime', frameDelay, imwriteOptions{:});
            else
                % Append to existing GIF
                imwrite(A, map, gifPath, 'gif', 'WriteMode', 'append', ...
                    'DelayTime', frameDelay, imwriteOptions{:});
            end
        end
        close(videoHandle);
    end
end
