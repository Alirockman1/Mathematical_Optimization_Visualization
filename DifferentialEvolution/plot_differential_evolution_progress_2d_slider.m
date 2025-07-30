function plot_differential_evolution_progress_2d_slider(optimizationData, varargin)
% PLOT_DIFFERENTIAL_EVOLUTION_PROGRESS_2D_SLIDER Visualizes 2D Differential Evolution optimization progress
%
%   PLOT_DIFFERENTIAL_EVOLUTION_PROGRESS_2D_SLIDER(OPTIMIZATIONDATA) creates
%   an interactive visualization of Differential Evolution optimization with:
%     - Four-panel view of DE operations (initialization, mutation, recombination, selection)
%     - Specialized 2D visualization with contour plots
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
%       'IncludeLegend' - Show legends (default: true)
%       'IncludeArrow' - Show mutation direction arrows (default: true)
%       'LegendOptions' - Legend customization options (cell array)
%       'ShowGrid' - Show grid lines ('on' or 'off', default: 'on')
%       Visual style parameters:
%           'MarkerEdgeColor', 'Marker', 'LineStyle', 'LineType', 
%           'MarkerFaceColor', 'MarkerSize', 'TitleFontSize', 'AxisFontSize', 'LegendFontSize'
%       'DesignVariableNames' - Custom names for design variables (cell array)
%       'UseDesignSpaceAsAxisLimits' - Use problem bounds for axes (default: true)
%       'PlotContour2DUseDesignSpaceLimit' - Use bounds for contour (default: true)
%       'PlotContour2DOptions' - Additional contour plot options (cell array)
%
%   OUTPUTS:
%       Interactive figure showing:
%       1. Four-panel view of DE operations:
%          - Initialization: Current population (blue)
%          - Mutation: base vectors (red) + differential vectors (arrows)
%          - Recombination: trial vectors (green) vs current population
%          - Selection: improved solutions (magenta) vs current population
%       2. Optional contour plots for 2D objective function visualization
%       3. Slider to navigate iterations
%       4. Optionally saves animation as GIF/video
%
%   KEY FEATURES:
%       - Optimized for 2D problems with clear visualization
%       - Comprehensive visualization of all DE operations
%       - Interactive exploration of optimization history
%       - Professional formatting with LaTeX labels
%       - Flexible saving options for publications
%       - Contour plot integration for objective function landscape
%
%   ALGORITHM DETAILS:
%       The visualization shows four key DE operations per iteration:
%       1. Initialization: Current population state
%       2. Mutation: base vectors + differential mutation vectors (with arrows)
%       3. Recombination: trial vectors vs current population elements
%       4. Selection: accepted solutions vs rejected solutions
%
%   REQUIREMENTS:
%       - Requires exactly 2 design variables for proper 2D visualization
%       - optimizationData must contain complete DE iteration history
%       - Objective function must be evaluable for contour generation
%
%   SEE ALSO:
%       plot_differential_evolution_progress_nd_slider, optimization_differential_evolution, 
%       plot_contour_objective_constraint_results_slider
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor, Contributor)
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
	parser.addParameter('SaveFolder',[]);
    parser.addParameter('SaveFigureOptions',{});
    parser.addParameter('CreateVideoGifIterations',true);
    parser.addParameter('FramesPerSecond',5);
    parser.addParameter('ImwriteOptions',{});
    parser.addParameter('VideoWriterProfile','Motion JPEG 2000');
    parser.addParameter('VideoWriterOptions',{'LosslessCompression',true});
    parser.addParameter('IncludeLegend',true);
    parser.addParameter('IncludeArrow',true);
    parser.addParameter('LegendOptions',{});
    parser.addParameter('ShowGrid', 'on', @(x) islogical(x) || ischar(x) && any(strcmpi(x, {'on', 'off'})));
    parser.addParameter('MarkerEdgeColor', 'b', @(x) ischar(x) || isstring(x));
    parser.addParameter('Marker', 'o', @(x) ischar(x) || isstring(x));
    parser.addParameter('LineStyle', '-', @(x) ischar(x) || isstring(x));
    parser.addParameter('LineType', 'line', @(x) ischar(x) || isstring(x));
    parser.addParameter('MarkerFaceColor', 'auto', @(x) ischar(x) || isstring(x));
    parser.addParameter('MarkerSize', 36, @(x) isnumeric(x) && isscalar(x));
    parser.addParameter('TitleFontSize', 14, @(x) isnumeric(x) && isscalar(x));
    parser.addParameter('AxisFontSize', 12, @(x) isnumeric(x) && isscalar(x));
    parser.addParameter('LegendFontSize', 12, @(x) isnumeric(x) && isscalar(x));
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
    handleCurrentPopulation = plot(ax1, NaN, NaN, '.', 'MarkerSize', 10);
    handleMutationInitialPoint = plot(ax2, NaN, NaN, '.', 'Color', 'r');
    handleMutationDirection = quiver(ax2, NaN, NaN, NaN,NaN,'AutoScale','on','Color','k');
    handleMutation1 = plot(ax2, NaN, NaN, 'o', 'Color', 'r');
    handleMutation2 = plot(ax3, NaN, NaN, '.', 'Color', 'r', 'MarkerSize', 10);
    handlePreviousPopulation1 = plot(ax3, NaN, NaN, '.', 'Color', 'b', 'MarkerSize', 10);
    handleRecombination2 = plot(ax4, NaN, NaN, '.', 'Color', 'g', 'MarkerSize', 10);
    handlePreviousPopulation2 = plot(ax4, NaN, NaN, '.', 'Color', 'b', 'MarkerSize', 10);

    % Titles, labels, limits for each subplot    for i = 1:4
        subplot(2, 2, i);
        title(subplotTitles{i});
        xlabel('x1'); ylabel('x2');
        xlim(xLimits);
        ylim(yLimits);
        grid(options.ShowGrid);

    if options.ShowContour
        % Create contour on a temporary axis (invisible)
        [x,y] = meshgrid(linspace(xLimits(1),xLimits(2),100),linspace(yLimits(1),yLimits(2),100));
        fullGrid = [x(:),y(:)];
        z = evaluate_optimization_objective(optimizationData.ProblemData.ObjectiveFunction,fullGrid);
        z = reshape(z, size(x));
        tempAxes = axes('Visible', 'off');
        [~, hContour] = contour(tempAxes, x, y, z, options.ContourLevels, options.ContourOptions{:});
        colormap('jet'); % Apply colormap globally if needed
        cbar = colorbar;
        cbar.Label.String = 'Objective Value';

        % Copy the contour to each of the 4 subplots
        copyobj(hContour, ax1);
        copyobj(hContour, ax2);
        copyobj(hContour, ax3);
        copyobj(hContour, ax4);

        % Delete the temporary axes
        delete(tempAxes);
    end

    % Create slider for iteration control
    slider = uicontrol('Style', 'slider', ...
        'Min', 1, 'Max', nIter, ...
        'Value', 1, 'SliderStep', [1/(nIter-1), 1/(nIter-1)], ...
        'Position', [100, 10, 400, 20], ...
        'Callback', @(src, event) updatePlot(round(src.Value)));

    % Function to update the plots based on slider position
    function updatePlot(i)
        % set(sliderText, 'String', sprintf('Iteration: %d', i));

        if(i==1)
			populationPrevious = optimizationData.InitialData.PopulationDesignInitial;
            updatePopulation = false(optimizationData.ProblemData.Options.PopulationSize,1);
		else
			populationPrevious = optimizationData.IterationData(i-1).PopulationDesignNew;
            updatePopulation = (optimizationData.IterationData(i).EvaluationTrialVectorObjectiveValue < ...
            optimizationData.IterationData(i-1).PopulationObjectiveValueNew);
        end
        mutatedPopulation = optimizationData.IterationData(i).MutationDonorVector;
        recombinedPopulation = optimizationData.IterationData(i).RecombinationTrialVector;
        selectedPopulation = optimizationData.IterationData(i).PopulationDesignNew;
        useDonorVector = optimizationData.IterationData(i).RecombinationUseDonorVector;

        sgtitle(sprintf('Algorithm: Differential Evolution, Population Size %d, Iteration %d', optimizationData.ProblemData.Options.PopulationSize, i)); % Adds overarching title

        % Current Population (Subplot 1)
        set(handleCurrentPopulation, 'XData', populationPrevious(:, 1), ...
                                     'YData', populationPrevious(:, 2));
        if options.IncludeLegend
            legend1 = legend(ax1, handleCurrentPopulation, {'Current Population'}, legendOptions{:});
            set(legend1, 'Box', 'on', 'Color', [0.9 0.9 0.9], 'EdgeColor', 'k');
        end

        % Mutation Step (Subplot 2)
        set(handleMutationInitialPoint, 'XData', optimizationData.IterationData(i).MutationInitialDesign(:,1), ...
                                        'YData', optimizationData.IterationData(i).MutationInitialDesign(:,2));
        if(options.IncludeArrow)
            set(handleMutationDirection, ...
            'XData', optimizationData.IterationData(i).MutationInitialDesign(:,1), ...
            'YData', optimizationData.IterationData(i).MutationInitialDesign(:,2), ...
            'UData', optimizationData.IterationData(i).MutationDirection(:,1), ...
            'VData', optimizationData.IterationData(i).MutationDirection(:,2),...
			'AutoScale','on','Color','k');
        end
        set(handleMutation1, 'XData', mutatedPopulation(:, 1), ...
                             'YData', mutatedPopulation(:, 2));
        if options.IncludeLegend
            legend2 = legend(ax2, [handleMutationInitialPoint, handleMutationDirection, handleMutation1], ...
                {'Mutation Initial Point', 'Mutation Direction', 'Mutant Population'}, legendOptions{:});
            set(legend2, 'Box', 'on', 'Color', [0.9 0.9 0.9], 'EdgeColor', 'k');
        end

        % Recombination Step (Subplot 3)
        set(handleMutation2, 'XData', mutatedPopulation(useDonorVector, 1), ...
                             'YData', mutatedPopulation(useDonorVector, 2));
        set(handlePreviousPopulation1, 'XData', populationPrevious(~useDonorVector, 1), ...
                                       'YData', populationPrevious(~useDonorVector, 2));
        if options.IncludeLegend
            legend3 = legend([handlePreviousPopulation1, handleMutation2], ...
                {'chosen from Current Population', 'chosen from Mutant Population'}, legendOptions{:});
            set(legend3, 'Box', 'on', 'Color', [0.9 0.9 0.9], 'EdgeColor', 'k');
        end

        % Selection Step (Subplot 4)
        set(handleRecombination2, 'XData', recombinedPopulation(updatePopulation, 1), ...
                                  'YData', recombinedPopulation(updatePopulation, 2));
        set(handlePreviousPopulation2, 'XData', populationPrevious(~updatePopulation, 1), ...
                                       'YData', populationPrevious(~updatePopulation, 2));
        if options.IncludeLegend
            legend4 = legend([handlePreviousPopulation2, handleRecombination2], ...
                {'chosen from Current Population', 'chosen from Recombined Population'}, legendOptions{:});
            set(legend4, 'Box', 'on', 'Color', [0.9 0.9 0.9], 'EdgeColor', 'k');
        end

        % Refresh the figure
        drawnow;
    end

    % Initialize the first iteration
    updatePlot(1);

    if options.CreateVideoGifIterations
        filename = [options.SaveFolder, 'IterationProgression'];
        videoHandle = VideoWriter(filename, options.VideoWriterProfile);
        for i = 1:2:length(videoWriterOptions)
            videoHandle.(videoWriterOptions{i}) = videoWriterOptions{i+1};
        end
        open(videoHandle);

        for i = 1:nIter
            slider.Value = i;
            updatePlot(i); % Update plot to match slider position
            currentFrame = getframe(figureHandle);
            writeVideo(videoHandle, currentFrame);
        end

        close(videoHandle);
    end

end
