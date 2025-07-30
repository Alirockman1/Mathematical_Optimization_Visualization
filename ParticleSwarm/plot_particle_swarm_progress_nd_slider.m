function plot_particle_swarm_progress_nd_slider(optimizationData, varargin)
% PLOT_PARTICLE_SWARM_PROGRESS_ND_SLIDER - Visualize the iterative progress 
% of a Particle Swarm Optimization (PSO) algorithm in n-dimensional design spaces 
% using 2D variable pair projections and a UI slider for interactive exploration.
%
% INPUTS:
%   optimizationData (struct)
%       Struct containing optimization run data, including problem definition,
%       design space bounds, initial population, velocities, and population 
%       positions over iterations. Expected fields include:
%           - InitialData.PopulationDesignInitial
%           - InitialData.OptimumCandidate
%           - IterationData(:).PopulationDesignNew
%           - IterationData(:).Velocities
%           - ProblemData.DesignSpaceLowerBound
%           - ProblemData.DesignSpaceUpperBound
%           - ProblemData.ObjectiveFunction
%           - ProblemData.Options.PopulationSize
%
%   varargin (Name-Value Pairs)
%       Optional arguments to customize the visualization:
%           'SaveFolder'                : Folder path to save outputs (default: [])
%           'SaveFigureOptions'         : Options for saving figures
%           'CreateVideoGifIterations'  : Generate video/gif of iterations (default: true)
%           'FramesPerSecond'           : FPS for video/gif output (default: 5)
%           'ImwriteOptions'            : Additional options for imwrite
%           'VideoWriterProfile'        : Profile for VideoWriter (default: 'Motion JPEG 2000')
%           'VideoWriterOptions'        : Additional VideoWriter properties
%           'ShowContour'               : Flag to show objective contours (default: true)
%           'ContourLevels'             : Number of contour levels (default: 10)
%           'ContourOptions'            : Additional contour plot options
%           'IncludeLegend'             : Include legends (default: true)
%           'IncludeArrow'              : Include velocity arrows (default: true)
%           'LegendOptions'             : Legend customizations
%           'ShowGrid'                  : Grid display ('on' or 'off') (default: 'on')
%           'MarkerEdgeColor'           : Edge color of particles (default: 'b')
%           'Marker'                    : Marker type for particles (default: 'o')
%           'LineStyle'                 : Line style for plot elements
%           'LineType'                  : Line type specifier
%           'MarkerFaceColor'           : Face color of markers
%           'MarkerSize'                : Size of particle markers (default: 36)
%           'AxisFontSize'              : Font size for axis labels (default: 12)
%           'LegendFontSize'            : Font size for legend (default: 12)
%           'DesignVariableNames'       : Cell array of design variable names (LaTeX strings)
%           'UseDesignSpaceAsAxisLimits': Use design bounds for axis limits (default: true)
%           'PlotContour2DUseDesignSpaceLimit': Flag for contour limit method (default: true)
%           'PlotContour2DOptions'      : Additional contour options for 2D plots
%
% OUTPUTS:
%   A MATLAB figure is displayed showing the evolution of the swarm across
%   iterations. If enabled, video and GIF outputs are saved in the specified folder.
%
% FUNCTIONALITY:
%   - Initializes a multi-panel plot where each subplot shows a 2D projection of the swarm.
%   - Updates plots with previous positions, velocities (as arrows), and new positions.
%   - Evaluates objective contours on 2D slices and overlays them if enabled.
%   - Includes a slider to navigate through optimization iterations interactively.
%   - Optionally saves iteration video and animated GIF of the progress.
%
% DEPENDENCIES:
%   - evaluate_optimization_objective
%   - merge_name_value_pair_argument
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor)
%   Copyright 2025 Ali Abbas Kapadia (Main Author) 
%   Copyright 2025 Gaurav Vaibhav (Contributor)
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

    % Problem info
    designVariables = size(optimizationData.InitialData.PopulationDesignInitial, 2);
    varPairs = nchoosek(1:designVariables, 2);
    nPairs = size(varPairs, 1);
    nIter = length(optimizationData.IterationData);

% --------------------------------------------------------------
% Define Plot layout
% --------------------------------------------------------------
    % Create figure
    figureHandle = figure('Name', 'Particle Swarm Progress', 'NumberTitle', 'off');
    set(gcf, 'Color', 'w');

    % Create subplots for all variable pairs in nPairs rows × 1 col
    axesHandles = gobjects(nPairs, 1);
    handlePreviousPopulation = gobjects(nPairs, 1);
    handleVelocities = gobjects(nPairs, 1);
    handleNewPopulation = gobjects(nPairs, 1);

    if designVariables == 2
        axesHandles = axes;

         % Set axis limits
        xLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(1), optimizationData.ProblemData.DesignSpaceUpperBound(1)];
        yLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(2), optimizationData.ProblemData.DesignSpaceUpperBound(2)];

        
        xlabel(axesHandles(1), options.DesignVariableNames{1}, 'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
        ylabel(axesHandles(1), options.DesignVariableNames{2}, 'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
        xlim(xLimits);
        ylim(yLimits);
        grid(options.ShowGrid);

        % Initialize plotting handles
        hold on;
        handlePreviousPopulation(1) = plot(NaN, NaN, '.', 'MarkerSize', 15, 'Color', 'b');
        handleVelocities(1) = quiver(NaN, NaN, NaN, NaN, 'AutoScale', 'on', 'Color', 'k');
        handleNewPopulation(1) = plot(NaN, NaN, 'o', 'Color', 'r');

        % Plot 2D contour slices
        if options.ShowContour
            % 1. Evaluate objective function grid
            [x, y] = meshgrid(linspace(xLimits(1), xLimits(2), 100), ...
                              linspace(yLimits(1), yLimits(2), 100));
            fullGrid = [x(:), y(:)];
            z = evaluate_optimization_objective(optimizationData.ProblemData.ObjectiveFunction, fullGrid);
            z = reshape(z, size(x));
        
            % 2. Create invisible figure + axis to generate contour objects
            [~, hContourLines] = contour(axesHandles, x, y, z, options.ContourLevels, options.ContourOptions{:});
            colormap('jet'); % Apply colormap globally if needed
            cbar = colorbar;
            cbar.Label.String = 'Objective Value';
        end

    elseif designVariables > 2
        % Subplot layout: max 3 columns
        numberColumns = 3;
        numberRows = ceil(nPairs / numberColumns);
        subplotMargin = 0.1;  % 10% spacing

        for pairIdx = 1:nPairs
            dimension1 = varPairs(pairIdx, 1);
            dimension2 = varPairs(pairIdx, 2);
            
            % Set axis limits
            xLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(dimension1), optimizationData.ProblemData.DesignSpaceUpperBound(dimension1)];
            yLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(dimension2), optimizationData.ProblemData.DesignSpaceUpperBound(dimension2)];
            
            row = floor((pairIdx-1)/numberColumns);
            column = mod((pairIdx-1), numberColumns);
    
            subplotWidth  = (1 - subplotMargin*(numberColumns+1)) / numberColumns;
            subplotHeight = (1 - subplotMargin*(numberRows+1)) / numberRows;
    
            left = subplotMargin + column * (subplotWidth + subplotMargin);
            bottom = 1 - subplotMargin - (row+1)*subplotHeight - row*subplotMargin;
    
            axesHandles(pairIdx) = axes('Position', [left, bottom, subplotWidth, subplotHeight]);
            hold(axesHandles(pairIdx), 'on');
    
            pair1 = varPairs(pairIdx, 1);
            pair2 = varPairs(pairIdx, 2);
    
            xlabel(axesHandles(pairIdx), options.DesignVariableNames{pair1}, 'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
            ylabel(axesHandles(pairIdx), options.DesignVariableNames{pair2}, 'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
            xlim(axesHandles(pairIdx), xLimits);
            ylim(axesHandles(pairIdx), yLimits);
            grid(axesHandles(pairIdx), options.ShowGrid);
    
            title(axesHandles(pairIdx), sprintf('%s vs %s', ...
                options.DesignVariableNames{pair2}, options.DesignVariableNames{pair1}), ...
                'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
    
            grid(axesHandles(pairIdx), options.ShowGrid);
            handlePreviousPopulation(pairIdx) = plot(axesHandles(pairIdx),NaN, NaN, '.', 'MarkerSize', 15,'Color','b');
            handleVelocities(pairIdx) = quiver(axesHandles(pairIdx),NaN, NaN, NaN, NaN, 'AutoScale','on','Color','k');
            handleNewPopulation(pairIdx) = plot(axesHandles(pairIdx),NaN, NaN, 'o', 'Color','r');

            % Plot 2D contour slices
            if options.ShowContour
                % 1. Evaluate objective function grid
                [x, y] = meshgrid(linspace(xLimits(1), xLimits(2), 100), ...
                                  linspace(yLimits(1), yLimits(2), 100));
                fullGrid = [x(:), y(:)];
                z = evaluate_optimization_objective(optimizationData.ProblemData.ObjectiveFunction, fullGrid);
                z = reshape(z, size(x));
            
                % 2. Create invisible figure + axis to generate contour objects
                [~, hContourLines] = contour(axesHandles(pairIdx), x, y, z, options.ContourLevels, options.ContourOptions{:});
                colormap('jet'); % Apply colormap globally if needed
                cbar = colorbar;
                cbar.Label.String = 'Objective Value';
            end
        end
    else
        error('At least two design variables are required for 2D plotting.');
    end

% --------------------------------------------------------------
% Define Slider UI element
% --------------------------------------------------------------
    if nIter > 1
        sliderStep = [1/(nIter-1), 1/(nIter-1)];
    else
        sliderStep = [1, 1]; % fallback for single iteration
    end

    slider = uicontrol('Style', 'slider', ...
        'Min', 1, 'Max', nIter, ...
        'Value', 1, ...
        'SliderStep', sliderStep, ...
        'Position', [100, 10, 400, 20], ...
        'Callback', @(src, ~) updatePlot(round(src.Value)));

% --------------------------------------------------------------
% Define Graph Update Function
% --------------------------------------------------------------
    function updatePlot(iter)
        % Get particle positions at iteration
        if(iter==1)
			populationPrevious = optimizationData.InitialData.PopulationDesignInitial;
		else
			populationPrevious = optimizationData.IterationData(iter-1).PopulationDesignNew;
        end

        % Update the title with iteration info
        sgtitle(sprintf('Particle Swarm Iteration %d, Population Size %d', iter, optimizationData.ProblemData.Options.PopulationSize)); % Adds overarching title

        for pairIdx = 1:nPairs
            pair1 = varPairs(pairIdx, 1);
            pair2 = varPairs(pairIdx, 2);
            set(handlePreviousPopulation(pairIdx), 'XData', populationPrevious(:, pair1), ...
                                     'YData', populationPrevious(:,pair2));
            set(handleVelocities(pairIdx), 'XData', populationPrevious(:,pair1),...
                                  'YData', populationPrevious(:,pair2),...
                                  'UData', optimizationData.IterationData(iter).Velocities(:,pair1),...
                                  'VData', optimizationData.IterationData(iter).Velocities(:,pair2));
            set(handleNewPopulation(pairIdx), 'XData', optimizationData.IterationData(iter).PopulationDesignNew(:,pair1),...
                                     'YData', optimizationData.IterationData(iter).PopulationDesignNew(:,pair2));
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

        for iter = 1:nIter
            slider.Value = iter;
            updatePlot(iter);

            pause(frameDelay);

            currentFrame = getframe(figureHandle);
            writeVideo(videoHandle, currentFrame);

            [A, map] = rgb2ind(frame2im(currentFrame), 256);
            if iter == 1
                imwrite(A, map, [filename, '.gif'], 'gif', 'LoopCount', inf, 'DelayTime', frameDelay, imwriteOptions{:});
            else
                imwrite(A, map, [filename, '.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', frameDelay, imwriteOptions{:});
            end

        end

        close(videoHandle);
    
    end
end