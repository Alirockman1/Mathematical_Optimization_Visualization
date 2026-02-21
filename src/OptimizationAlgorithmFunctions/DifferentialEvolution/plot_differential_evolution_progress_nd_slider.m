function plot_differential_evolution_progress_nd_slider(optimizationData, varargin)
% PLOT_DIFFERENTIAL_EVOLUTION_PROGRESS_ND_SLIDER Visualizes Differential Evolution optimization progress
%
%   PLOT_DIFFERENTIAL_EVOLUTION_PROGRESS_ND_SLIDER(OPTIMIZATIONDATA) creates
%   an interactive visualization of Differential Evolution optimization with:
%     - Multi-panel view of DE operations (mutation, recombination, selection)
%     - 2D projections for high-dimensional problems
%     - Interactive slider to navigate iterations
%     - Option to save as GIF/video
%
%   INPUTS:
%       optimizationData - Structure containing DE optimization results with fields:
%           .ProblemData - Problem definition including bounds and objective
%           .IterationData - Array of structs with optimization history
%           .InitialData - Initial population and values
%
%   OPTIONAL PARAMETERS (name-value pairs):
%       'SaveFolder' - Path to save output figures/videos (default: [])
%       'SaveFigureOptions' - Options for figure saving (cell array)
%       'CreateVideoGifIterations' - Enable GIF/video creation (default: true)
%       'FramesPerSecond' - Animation frame rate (default: 5)
%       'ImwriteOptions' - Options for GIF creation (cell array)
%       'VideoWriterProfile' - Video profile for VideoWriter (default: 'Motion JPEG 2000')
%       'VideoWriterOptions' - Video writing options (cell array)
%       'ShowContour' - Show objective function contours (default: true)
%       'ContourLevels' - Number of contour levels (default: 10)
%       'ContourOptions' - Additional contour plot options (cell array)
%       'IncludeLegend' - Show legends (default: true)
%       'IncludeArrow' - Show mutation direction arrows (default: true)
%       'LegendOptions' - Legend customization options (cell array)
%       'ShowGrid' - Show grid lines ('on' or 'off', default: 'on')
%       Visual style parameters:
%           'MarkerEdgeColor', 'Marker', 'LineStyle', 'LineType', 
%           'MarkerFaceColor', 'MarkerSize', 'AxisFontSize', 'LegendFontSize'
%       'DesignVariableNames' - Custom names for design variables (cell array)
%       'UseDesignSpaceAsAxisLimits' - Use problem bounds for axes (default: true)
%       'PlotContour2DUseDesignSpaceLimit' - Use bounds for contour (default: true)
%       'PlotContour2DOptions' - Additional contour plot options (cell array)
%
%   OUTPUTS:
%       Interactive figure showing:
%       1. Four-panel view of DE operations:
%          - Current population
%          - Mutation step with direction vectors
%          - Recombination step
%          - Selection step
%       2. For ND problems: grid of 2D projections
%       3. Slider to navigate iterations
%       4. Optionally saves animation as GIF/video
%
%   KEY FEATURES:
%       - Comprehensive visualization of all DE operations
%       - Handles both 2D and high-dimensional problems
%       - Interactive exploration of optimization history
%       - Professional formatting with LaTeX labels
%       - Flexible saving options for publications
%
%   ALGORITHM DETAILS:
%       The visualization shows four key DE operations per iteration:
%       1. Current population (blue)
%       2. Mutation: base vectors (red) + differential vectors (arrows)
%       3. Recombination: trial vectors (green) vs current population
%       4. Selection: improved solutions (magenta) vs current population
%
%   SEE ALSO:
%       differential_evolution_optimization, plot_contour_objective_constraint_results_slider
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor)
%   Copyright 2025 Ali Abbas Kapadia (Main Author) 
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

    % Define the color sequences for each DE step (1 to 4)
    colorSequences = {
        {'b'},             % Step 1: Current Population (initial/previous) - Blue
        {'b', 'g'},        % Step 2: Mutation (prev + mutant) - Blue + Green
        {'g', 'k'},        % Step 3: Recombination (mutant + trial) - Green + Black
        {'b', 'm'}         % Step 4: Selection (prev + selected) - Blue + Magenta
    };

    % Problem info
    designVariables = size(optimizationData.InitialData.PopulationDesignInitial, 2);
    varPairs = nchoosek(1:designVariables, 2);
    nPairs = size(varPairs, 1);
    nIter = length(optimizationData.IterationData);
    subplotTitles = {'Initialization  Step', 'Mutation Step', 'Recombination Step', 'Selection Step'};

    % --------------------------------------------------------------
    % Define Plot layout
    % --------------------------------------------------------------
    % Create figure
    figureHandle = figure('Name', 'Differential evolution progress', 'NumberTitle', 'off');
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

        handleCurrentPopulation = plot(ax1, NaN, NaN, '.', 'MarkerSize', 10);
        handleMutationInitialPoint = plot(ax2, NaN, NaN, '.', 'Color', 'r');
        handleMutationDirection = quiver(ax2, NaN, NaN, NaN,NaN,'AutoScale','on','Color','k');
        handleMutation1 = plot(ax2, NaN, NaN, 'o', 'Color', 'r');
        handleMutation2 = plot(ax3, NaN, NaN, '.', 'Color', 'r', 'MarkerSize', 10);
        handlePreviousPopulation1 = plot(ax3, NaN, NaN, '.', 'Color', 'b', 'MarkerSize', 10);
        handleRecombination2 = plot(ax4, NaN, NaN, '.', 'Color', 'g', 'MarkerSize', 10);
        handlePreviousPopulation2 = plot(ax4, NaN, NaN, '.', 'Color', 'b', 'MarkerSize', 10);

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

            handleCurrentPopulation(pairIdx)   = plot(axesHandles(pairIdx,1), NaN, NaN, '.', 'MarkerSize', 10);
            handleMutationInitialPoint(pairIdx) = plot(axesHandles(pairIdx,2), NaN, NaN, '.', 'Color', 'r');
            handleMutationDirection(pairIdx)   = quiver(axesHandles(pairIdx,2), NaN, NaN, NaN,NaN,'AutoScale','on','Color','k');
            handleMutation1(pairIdx)           = plot(axesHandles(pairIdx,2), NaN, NaN, 'o', 'Color', 'r');
            handleMutation2(pairIdx)           = plot(axesHandles(pairIdx,3), NaN, NaN, '.', 'Color', 'r', 'MarkerSize', 10);
            handlePreviousPopulation1(pairIdx) = plot(axesHandles(pairIdx,3), NaN, NaN, '.', 'Color', 'b', 'MarkerSize', 10);
            handleRecombination2(pairIdx)      = plot(axesHandles(pairIdx,4), NaN, NaN, '.', 'Color', 'g', 'MarkerSize', 10);
            handlePreviousPopulation2(pairIdx) = plot(axesHandles(pairIdx,4), NaN, NaN, '.', 'Color', 'b', 'MarkerSize', 10);

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

% --------------------------------------------------------------
% Define Slider UI element
% --------------------------------------------------------------
    % Create slider for iteration control
    slider = uicontrol('Style', 'slider', ...
        'Min', 1, 'Max', nIter, ...
        'Value', 1, 'SliderStep', [1/(nIter-1), 1/(nIter-1)], ...
        'Position', [100, 10, 400, 20], ...
        'Callback', @(src, event) updatePlot(round(src.Value)));

    % Function to update the plots based on slider position
    function updatePlot(iter)
        if(iter==1)
			populationPrevious = optimizationData.InitialData.PopulationDesignInitial;
            updatePopulation = false(optimizationData.ProblemData.Options.PopulationSize,1);
		else
			populationPrevious = optimizationData.IterationData(iter-1).PopulationDesignNew;
            updatePopulation = (optimizationData.IterationData(iter).EvaluationTrialVectorObjectiveValue < ...
            optimizationData.IterationData(iter-1).PopulationObjectiveValueNew);
        end
        mutatedPopulation = optimizationData.IterationData(iter).MutationDonorVector;
        recombinedPopulation = optimizationData.IterationData(iter).RecombinationTrialVector;
        selectedPopulation = optimizationData.IterationData(iter).PopulationDesignNew;
        useDonorVector = optimizationData.IterationData(iter).RecombinationUseDonorVector;

        sgtitle(sprintf('Differential Evolution Iteration %d, Population Size %d',  iter, optimizationData.ProblemData.Options.PopulationSize)); % Adds overarching title
        
        for plotIter = 1:nPairs
            i1 = varPairs(plotIter, 1);
            i2 = varPairs(plotIter, 2);
    
            % Current Population (Subplot 1)
            ax1 = axesHandles(plotIter, 1);
            set(handleCurrentPopulation(plotIter), 'XData', populationPrevious(:, i1), ...
                                        'YData', populationPrevious(:, i2));
            if options.IncludeLegend
                legend1 = legend(ax1, handleCurrentPopulation(plotIter), {'Current Population'}, legendOptions{:});
                set(legend1, 'Box', 'on', 'Color', [0.9 0.9 0.9], 'EdgeColor', 'k');
            end

            % Mutation Step (Subplot 2)\
            ax2 = axesHandles(plotIter, 2);
            set(handleMutationInitialPoint(plotIter), 'XData', optimizationData.IterationData(iter).MutationInitialDesign(:,i1), ...
                                            'YData', optimizationData.IterationData(iter).MutationInitialDesign(:,i2));
            if(options.IncludeArrow)
                set(handleMutationDirection(plotIter), ...
                'XData', optimizationData.IterationData(iter).MutationInitialDesign(:,i1), ...
                'YData', optimizationData.IterationData(iter).MutationInitialDesign(:,i2), ...
                'UData', optimizationData.IterationData(iter).MutationDirection(:,i1), ...
                'VData', optimizationData.IterationData(iter).MutationDirection(:,i2),...
			    'AutoScale','on','Color','k');
            end
            set(handleMutation1(plotIter), 'XData', mutatedPopulation(:, i1), ...
                                 'YData', mutatedPopulation(:, i2));
            if options.IncludeLegend
                legend2 = legend(ax2, [handleMutationInitialPoint(plotIter), handleMutationDirection(plotIter), handleMutation1(plotIter)], ...
                    {'Mutation Initial Point', 'Mutation Direction', 'Mutant Population'}, legendOptions{:});
                set(legend2, 'Box', 'on', 'Color', [0.9 0.9 0.9], 'EdgeColor', 'k');
            end
    
            % Recombination Step (Subplot 3)
            ax3 = axesHandles(plotIter, 3);
            set(handleMutation2(plotIter), 'XData', mutatedPopulation(useDonorVector, i1), ...
                                 'YData', mutatedPopulation(useDonorVector, i2));
            set(handlePreviousPopulation1(plotIter), 'XData', populationPrevious(~useDonorVector, i1), ...
                                           'YData', populationPrevious(~useDonorVector, i2));
            if options.IncludeLegend
                legend3 = legend([handlePreviousPopulation1(plotIter), handleMutation2(plotIter)], ...
                    {'chosen from Current Population', 'chosen from Mutant Population'}, legendOptions{:});
                set(legend3, 'Box', 'on', 'Color', [0.9 0.9 0.9], 'EdgeColor', 'k');
            end
    
            % Selection Step (Subplot 4)
            ax4 = axesHandles(plotIter, 4);
            set(handleRecombination2(plotIter), 'XData', recombinedPopulation(updatePopulation, i1), ...
                                      'YData', recombinedPopulation(updatePopulation, i2));
            set(handlePreviousPopulation2(plotIter), 'XData', populationPrevious(~updatePopulation, i1), ...
                                           'YData', populationPrevious(~updatePopulation, i2));
            if options.IncludeLegend
                legend4 = legend([handlePreviousPopulation2(plotIter), handleRecombination2(plotIter)], ...
                    {'chosen from Current Population', 'chosen from Recombined Population'}, legendOptions{:});
                set(legend4, 'Box', 'on', 'Color', [0.9 0.9 0.9], 'EdgeColor', 'k');
            end
        end
        % Refresh the figure
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
        for i = 1:2:length(videoWriterOptions)
            videoHandle.(videoWriterOptions{i}) = videoWriterOptions{i+1};
        end

        open(videoHandle);

        for playbackiter = 1:nIter
            slider.Value = playbackiter;
            updatePlot(playbackiter); % Update plot to match slider position
            
            pause(frameDelay);

            currentFrame = getframe(figureHandle);
            writeVideo(videoHandle, currentFrame);

            [A, map] = rgb2ind(frame2im(currentFrame), 256);
            if playbackiter == 1
                imwrite(A, map, [filename, '.gif'], 'gif', 'LoopCount',...
                    inf, 'DelayTime', frameDelay, imwriteOptions{:});
            else
                imwrite(A, map, [filename, '.gif'], 'gif', 'WriteMode',...
                    'append', 'DelayTime', frameDelay, imwriteOptions{:});
            end
        end

        close(videoHandle);
    end
end
