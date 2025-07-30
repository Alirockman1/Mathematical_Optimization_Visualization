function plot_genetic_algorithm_progress_2d_slider(optimizationData, varargin)
%PLOT_GENETIC_ALGORITHM_PROGRESS_2D_SLIDER Visualize GA progress in 2D design space with interactive slider.
%
%   PLOT_GENETIC_ALGORITHM_PROGRESS_2D_SLIDER(OPTIMIZATIONDATA, OPTIONS) 
%       creates an interactive visualization of genetic algorithm progress 
%       in 2D design space, showing the evolution of different GA steps 
%       (Current Population, Selection, Crossover, Mutation) across 
%       iterations using four subplots and a slider for iteration control.
%
%   INPUTS:
%       OPTIMIZATIONDATA    - Struct containing GA optimization results:
%           .ProblemData    - Problem definition with design space bounds
%           .InitialData    - Initial population data
%           .IterationData  - Array of iteration data containing:
%               .SelectionIndex      - Indices of selected individuals
%               .MutationIsMutated   - Boolean array indicating mutations
%               .SelectionPopulation - Selected parent population
%               .CrossoverPopulation - Crossover offspring population
%               .MutationPopulation  - Final mutated population
%
%   NAME-VALUE PAIR ARGUMENTS:
%       'SaveFolder'            - String path for saving output files (default: [])
%       'SaveFigureOptions'     - Cell array of figure save options (default: {})
%       'CreateVideoGifIterations' - Logical flag to create video/GIF (default: true)
%       'FramesPerSecond'       - Video frame rate (default: 5)
%       'ImwriteOptions'        - Cell array of image write options (default: {})
%       'VideoWriterProfile'    - Video codec profile (default: 'Motion JPEG 2000')
%       'VideoWriterOptions'    - Cell array of video writer options (default: {'LosslessCompression',true})
%       'ShowContour'           - Logical flag to show objective contours (default: true)
%       'ContourLevels'         - Number of contour levels (default: 10)
%       'ContourOptions'        - Cell array of contour plot options (default: {})
%       'IncludeLegend'         - Logical flag to include legends (default: true)
%       'IncludeArrow'          - Logical flag to include arrows (default: true)
%       'LegendOptions'         - Cell array of legend options (default: {})
%       'ShowGrid'              - Grid display option 'on'/'off' (default: 'on')
%       'MarkerEdgeColor'       - Marker edge color (default: 'b')
%       'Marker'                - Marker style (default: 'o')
%       'LineStyle'             - Line style (default: '-')
%       'LineType'              - Line type (default: 'line')
%       'MarkerFaceColor'       - Marker face color (default: 'auto')
%       'MarkerSize'            - Marker size (default: 36)
%       'AxisFontSize'          - Font size for axis labels (default: 12)
%       'LegendFontSize'        - Font size for legends (default: 12)
%
%   OUTPUTS:
%       None. Interactive figure is displayed with four subplots showing:
%       1. Current Population - Previous generation individuals
%       2. Selection Step - Current population with selected parents highlighted
%       3. Crossover Step - Selected parents with crossover offspring
%       4. Mutation Step - Crossover population with mutated individuals
%
%   FUNCTIONALITY:
%       1. Creates 2x2 subplot layout for GA visualization steps
%       2. Displays objective function contours as background (optional)
%       3. Uses different colors to distinguish GA populations at each step
%       4. Provides interactive slider to navigate through iterations
%       5. Supports video/GIF creation for animation export
%       6. Fully customizable visual appearance via options
%       7. Real-time plot updates with efficient handle graphics
%
%   NOTES:
%       - Requires 2D design space (exactly 2 design variables)
%       - Contour visualization requires objective function evaluation
%       - Video creation saves both .avi and .gif formats
%       - Optimized for educational demonstration of GA mechanics
%       - Compatible with standard GA optimization data structures
%
%   EXAMPLE:
%       % After running genetic algorithm optimization
%       plot_genetic_algorithm_progress_2d_slider(gaResults, ...
%           'ShowContour', true, ...
%           'ContourLevels', 15, ...
%           'CreateVideoGifIterations', false);
%
%   SEE ALSO:
%       plot_genetic_algorithm_progress_nd_slider, optimization_genetic_algorithm,
%       plot_objective_value_progress_2d
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

    % Parse input parameters
    parser = inputParser;
	parser.addParameter('SaveFolder',[]);
    parser.addParameter('SaveFigureOptions',{});
    parser.addParameter('CreateVideoGifIterations',true);
    parser.addParameter('FramesPerSecond',5);
    parser.addParameter('ImwriteOptions',{});
    parser.addParameter('VideoWriterProfile','Motion JPEG 2000');
    parser.addParameter('VideoWriterOptions',{'LosslessCompression',true});
    parser.addParameter('ShowContour', true);
    parser.addParameter('ContourLevels',10);
    parser.addParameter('ContourOptions',{});
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
    parser.addParameter('AxisFontSize', 12, @(x) isnumeric(x) && isscalar(x));
    parser.addParameter('LegendFontSize', 12, @(x) isnumeric(x) && isscalar(x));   
    parser.parse(varargin{:});
    options = parser.Results;

    defaultVideoWriterOptions = {'FrameRate',options.FramesPerSecond};
    [~,videoWriterOptions] = merge_name_value_pair_argument(defaultVideoWriterOptions,options.VideoWriterOptions);
    defaultImwriteOptions = {};
    [~,imwriteOptions] = merge_name_value_pair_argument(defaultImwriteOptions,options.ImwriteOptions);
    frameDelay = 1/options.FramesPerSecond;
    defaultLegendOptions = {'location','southeast','FontSize', 6};
    [~,legendOptions] = merge_name_value_pair_argument(defaultLegendOptions,options.LegendOptions);

    % Number of iterations
    nIter = length(optimizationData.IterationData);

    % Create figure and axis
    figureHandle = figure('Name', 'Genetic Algorithm Progress', 'NumberTitle', 'off');

    % Set axis limits
    xLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(1), optimizationData.ProblemData.DesignSpaceUpperBound(1)];
    yLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(2), optimizationData.ProblemData.DesignSpaceUpperBound(2)];

    % Generate meshgrid for contours
    if options.ShowContour
        [x,y] = meshgrid(linspace(xLimits(1),xLimits(2),100),linspace(yLimits(1),yLimits(2),100));
        fullGrid = [x(:),y(:)];
        z = evaluate_optimization_objective(optimizationData.ProblemData.ObjectiveFunction,fullGrid);
        z = reshape(z, size(x));
        % Set colormap and global color axis
        colormap('jet');
        cmin = min(z(:));
        cmax = max(z(:));
        clim([cmin cmax]); % Same color range for all subplots
    end

    % Create subplots for the 4 steps
    ax1 = subplot(2, 2, 1); hold on;
    ax2 = subplot(2, 2, 2); hold on;
    ax3 = subplot(2, 2, 3); hold on;
    ax4 = subplot(2, 2, 4); hold on;
    set(gcf, 'Color', 'w');  

    % Create empty plot handles to update later
    handleCurrentPopulation1 = plot(ax1, NaN, NaN, '.', 'MarkerSize', 10,'Color','b');
    handleCurrentPopulation2 = plot(ax2, NaN, NaN, '.', 'MarkerSize', 10,'Color','b');
    handleSelection1 = plot(ax2, NaN, NaN, '.', 'MarkerSize', 10, 'Color', 'r');
    handleSelection2 = plot(ax3, NaN, NaN, '.', 'MarkerSize', 10, 'Color', 'r');
    handleCrossOver1 = plot(ax3, NaN, NaN, '.', 'MarkerSize', 10, 'Color', 'k');
    handleCrossOver2 = plot(ax4, NaN, NaN, '.', 'MarkerSize', 10, 'Color', 'k');
    handleMutation = plot(ax4, NaN, NaN, '.', 'MarkerSize', 10, 'Color', 'g');

    % Titles, labels, limits for each subplot
    subplotTitles = {'Current Population', 'Selection Step', 'CrossOver Step', 'Mutation Step'};
    for i = 1:4
        subplot(2, 2, i);
        title(subplotTitles{i});
        xlabel('x1','interpreter','latex','FontSize',options.AxisFontSize);
        ylabel('x2','interpreter','latex','FontSize',options.AxisFontSize);
        xlim(xLimits);
        ylim(yLimits);
        grid(options.ShowGrid);
    end

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
        else
            populationPrevious = optimizationData.IterationData(i-1).MutationPopulation;
        end
        selectionIndex = optimizationData.IterationData(i).SelectionIndex;
        isMutate = optimizationData.IterationData(i).MutationIsMutated;
        selectedPopulation = optimizationData.IterationData(i).SelectionPopulation;
        crossOverPopulation = optimizationData.IterationData(i).CrossoverPopulation;
        mutatedPopulation = optimizationData.IterationData(i).MutationPopulation;

        sgtitle(sprintf('Algorithm: Genetic Algo Iteration %d, Population Size %d', i, optimizationData.ProblemData.Options.PopulationSize)); % Adds overarching title

        % Current Population (Subplot 1)
        set(handleCurrentPopulation1, 'XData', populationPrevious(:, 1), ...
                                     'YData', populationPrevious(:, 2));
        if options.IncludeLegend
            legend1 = legend(ax1, handleCurrentPopulation1, {'Current Population'}, legendOptions{:});
            set(legend1, 'Box', 'on', 'Color', [0.9 0.9 0.9], 'EdgeColor', 'k');
        end

        % Mutation Step (Subplot 2)
        set(handleCurrentPopulation2, 'XData', populationPrevious(:, 1), ...
                                     'YData', populationPrevious(:, 2));
        set(handleSelection1, 'XData', selectedPopulation(selectionIndex,1), ...
                             'YData', selectedPopulation(selectionIndex,2));

        if options.IncludeLegend
            legend2 = legend(ax2, [handleCurrentPopulation2, handleSelection1], ...
                {'Current Pop', 'Selected from Current Pop'}, legendOptions{:});
            set(legend2, 'Box', 'on', 'Color', [0.9 0.9 0.9], 'EdgeColor', 'k');
        end

        % Recombination Step (Subplot 3)
        set(handleSelection2, 'XData', selectedPopulation(:,1), ...
                             'YData', selectedPopulation(:,2));
        set(handleCrossOver1, 'XData', crossOverPopulation(:, 1), ...
                              'YData', crossOverPopulation(:, 2));
        if options.IncludeLegend
            legend3 = legend([handleSelection2, handleCrossOver1], ...
                {'Selected Population', 'CrossOver Population'}, legendOptions{:});
            set(legend3, 'Box', 'on', 'Color', [0.9 0.9 0.9], 'EdgeColor', 'k');
        end

        % Selection Step (Subplot 4)
        set(handleCrossOver2, 'XData', crossOverPopulation(~isMutate,1), ...
                              'YData', crossOverPopulation(~isMutate,2));
        set(handleMutation, 'XData', mutatedPopulation(isMutate,1), ...
                            'YData', mutatedPopulation(isMutate,2));
        if options.IncludeLegend
            legend4 = legend([handleCrossOver2, handleMutation], ...
                {'chosen from CrossOver Pop', 'Mutated Population'}, legendOptions{:});
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


