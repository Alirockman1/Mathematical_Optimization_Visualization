function plot_particle_swarm_progress_2d_slider(optimizationData, varargin)
% PLOT_PARTICLE_SWARM_PROGRESS_2D_SLIDER - Visualize the iterative progress 
% of a Particle Swarm Optimization (PSO) algorithm in 2-dimensional design spaces 
% using a UI slider for interactive exploration.
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
%           'SaveFolder'                        : Folder path to save outputs (default: [])
%           'SaveFigureOptions'                 : Options for saving figures
%           'CreateVideoGifIterations'          : Generate video/gif of iterations (default: true)
%           'FramesPerSecond'                   : FPS for video/gif output (default: 5)
%           'ImwriteOptions'                    : Additional options for imwrite
%           'VideoWriterProfile'                : Profile for VideoWriter (default: 'Motion JPEG 2000')
%           'VideoWriterOptions'                : Additional VideoWriter properties
%           'IncludeLegend'                     : Include legends (default: true)
%           'IncludeArrow'                      : Include velocity arrows (default: true)
%           'LegendOptions'                     : Legend customizations
%           'ShowGrid'                          : Grid display ('on' or 'off') (default: 'on')
%           'MarkerEdgeColor'                   : Edge color of particles (default: 'b')
%           'Marker'                            : Marker type for particles (default: 'o')
%           'LineStyle'                         : Line style for plot elements (default: '-')
%           'LineType'                          : Line type specifier (default: 'line')
%           'MarkerFaceColor'                   : Face color of markers (default: 'auto')
%           'MarkerSize'                        : Size of particle markers (default: 36)
%           'TitleFontSize'                     : Font size for titles (default: 14)
%           'AxisFontSize'                      : Font size for axis labels (default: 12)
%           'LegendFontSize'                    : Font size for legend (default: 12)
%           'DesignVariableNames'               : Cell array of design variable names (LaTeX strings)
%           'UseDesignSpaceAsAxisLimits'        : Use design bounds for axis limits (default: true)
%           'PlotContour2DUseDesignSpaceLimit'  : Flag for contour limit method (default: true)
%           'PlotContour2DOptions'              : Additional contour options for 2D plots
%
% OUTPUTS:
%   A MATLAB figure is displayed showing the evolution of the swarm across
%   iterations in 2D space. If enabled, video and GIF outputs are saved in the specified folder.
%
% FUNCTIONALITY:
%   - Initializes a 2D plot showing the particle swarm evolution over iterations.
%   - Updates plots with previous positions, velocities (as arrows), and new positions.
%   - Evaluates and overlays objective function contours for visualization context.
%   - Includes a slider to navigate through optimization iterations interactively.
%   - Optionally saves iteration video and animated GIF of the progress.
%   - Displays legend with color-coded particle states and velocity vectors.
%
% DEPENDENCIES:
%   - plot_point_optimization_contour_2d
%   - merge_name_value_pair_argument
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor)
%   Copyright 2025 Ali Abbas Kapadia (Main Author) 
%   Copyright 2025 Gaurav Vaibhav (Contributor)
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
    figureHandle = figure('Name', 'Particle swarm Progress', 'NumberTitle', 'off');
    ax = axes;
    set(gcf, 'Color', 'w');  
    hold(ax,'on');

    % Set axis limits
    xLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(1), optimizationData.ProblemData.DesignSpaceUpperBound(1)];
    yLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(2), optimizationData.ProblemData.DesignSpaceUpperBound(2)];
    xlim(xLimits);
    ylim(yLimits);

    % Create empty plot handles to update later
    handlePreviousPopulation = plot(NaN, NaN, '.', 'MarkerSize', 15,'Color','b');
    handleVelocities = quiver(NaN, NaN, NaN, NaN, 'AutoScale','on','Color','k');
    handleNewPopulation = plot(NaN, NaN, 'o', 'Color','r');

     % Plot contour (2D only)
    if length(optimizationData.InitialData.OptimumCandidate) == 2
        if options.PlotContour2DUseDesignSpaceLimit
            contourDesignSpace = [optimizationData.ProblemData(1).DesignSpaceLowerBound;
                                  optimizationData.ProblemData(1).DesignSpaceUpperBound];
        else
            optimumCandidate = cell2mat(arrayfun(@(d) d.OptimumCandidate, optimizationData.IterationData, 'UniformOutput', false));
            contourDesignSpace = [min(optimumCandidate,[],1); max(optimumCandidate,[],1)];
        end
        hContour = plot_point_optimization_contour_2d(figureHandle, optimizationData.ProblemData(1).ObjectiveFunction, contourDesignSpace, ...
            'DesignVariableNames', options.DesignVariableNames, ...
            'AxisFontSize', options.AxisFontSize, ...
            options.PlotContour2DOptions{:});
        if isgraphics(hContour)
            set(findall(hContour, '-property', 'HandleVisibility'), 'HandleVisibility', 'off');
        end
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
			populationPrevious = optimizationData.IterationData(i-1).PopulationDesignNew;
        end
       
        sgtitle(sprintf('Algorithm: Particle Swarm Iteration %d, Population Size %d', i, optimizationData.ProblemData.Options.PopulationSize)); % Adds overarching title

        set(handlePreviousPopulation, 'XData', populationPrevious(:, 1), ...
                                     'YData', populationPrevious(:, 2));
        set(handleVelocities, 'XData', populationPrevious(:,1),...
                              'YData', populationPrevious(:,2),...
                              'UData', optimizationData.IterationData(i).Velocities(:,1),...
                              'VData', optimizationData.IterationData(i).Velocities(:,2));
        set(handleNewPopulation, 'XData', optimizationData.IterationData(i).PopulationDesignNew(:,1),...
                                 'YData', optimizationData.IterationData(i).PopulationDesignNew(:,2));

        if options.IncludeLegend
            legend1 = legend(ax, [handlePreviousPopulation, handleVelocities, handleNewPopulation],...
                {'Old Population','Population Velocity', 'New Population'}, legendOptions{:});
            set(legend1, 'Box', 'on', 'Color', [0.9 0.9 0.9], 'EdgeColor', 'k');
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