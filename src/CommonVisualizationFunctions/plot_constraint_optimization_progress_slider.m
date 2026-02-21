function plot_constraint_optimization_progress_slider(optimizationData, objectiveFunction, varargin)
%PLOT_CONSTRAINT_OPTIMIZATION_PROGRESS_SLIDER Visualizes constrained optimization progress.
%
%   PLOT_CONSTRAINT_OPTIMIZATION_PROGRESS_SLIDER(OPTIMIZATIONDATA, OBJECTIVEFUNCTION, ...)
%   creates an interactive visualization of constrained optimization progress with:
%       - Contour plot of objective function (2D only)
%       - Iteration path animation
%       - Interactive slider control
%       - Optional GIF/video export
%
%   INPUTS:
%       OPTIMIZATIONDATA - Structure containing optimization history:
%           .ProblemData - Problem definition with bounds
%           .InitialData - Initial point information
%           .IterationData - Array of iteration records
%       OBJECTIVEFUNCTION - Function handle for objective evaluation
%
%   NAME-VALUE PAIR ARGUMENTS:
%       'OptimizationMethod' - Name of method used (default: 'Conjugate gradient')
%       'SaveFolder' - Directory to save outputs (default: [])
%       'CreateVideoGifIterations' - Enable animation export (default: true)
%       'FramesPerSecond' - Animation frame rate (default: 5)
%       'DesignVariableNames' - Custom names for variables (default: x₁,x₂,...)
%       'ArrowScale' - Scaling factor for direction arrows (default: 0.01)
%       (Complete list of 15+ parameters available in function)
%
%   FUNCTIONALITY:
%       1. Generates contour plot of objective function (2D problems)
%       2. Plots optimization path with iteration markers
%       3. Provides interactive slider to explore iterations
%       4. Supports animation export as GIF/video
%       5. Maintains design space aspect ratio
%
%   NOTES:
%       - Designed for 2D optimization problems
%       - Uses Arial font consistently throughout
%       - Automatically handles design space bounds
%       - Creates high-quality output for presentations/publications
%
%   EXAMPLE:
%       % Visualize optimization progress
%       plot_constraint_optimization_progress_slider(optData, @objFunc, ...
%           'SaveFolder', 'results/', 'FramesPerSecond', 10)
%
%   SEE ALSO:
%       plot_point_optimization_contour_2d, optimization_active_set
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor)
%   Copyright 2025 Ali Abbas Kapadia (Author) 
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

    parser = inputParser;
    parser.addParameter('OptimizationMethod', 'Conjugate gradient', @(x) ischar(x) || isstring(x));
    parser.addParameter('SaveFolder', []);
    parser.addParameter('SaveFigureOptions', {});
    parser.addParameter('CreateVideoGifIterations', true);
    parser.addParameter('FramesPerSecond', 5);
    parser.addParameter('ImwriteOptions', {});
    parser.addParameter('VideoWriterProfile', 'Motion JPEG 2000');
    parser.addParameter('VideoWriterOptions', {'LosslessCompression', true});
    parser.addParameter('TitleFontSize', 14, @(x) isnumeric(x) && isscalar(x));
    parser.addParameter('AxisFontSize', 12, @(x) isnumeric(x) && isscalar(x));
    parser.addParameter('DesignVariableNames',{});
    parser.addParameter('UseDesignSpaceAsAxisLimits',true);
    parser.addParameter('PlotContour2DUseDesignSpaceLimit',true);
    parser.addParameter('PlotContour2DOptions',{});
    parser.addParameter('ArrowScale',0.01, @(x) isnumeric(x) && isscalar(x));
    parser.parse(varargin{:});
    options = parser.Results;

    set(groot, ...
    'DefaultAxesFontName', 'Arial', ...
    'DefaultAxesFontSize', options.AxisFontSize);

    % Setup design variable names
    if isempty(options.DesignVariableNames)
        for i = 1:length(optimizationData.InitialData.OptimumCandidate)
            options.DesignVariableNames{i} = sprintf('$$x_{%d}$$', i);
        end
    end

    % Setup video/gif saving
    frameDelay = 1 / options.FramesPerSecond;
    gif_filename = fullfile(options.SaveFolder, 'Optimization_progress.gif');

    % Create figure
    figureHandle = figure;
    set(gcf, 'Color', 'w');
    set(gca, 'FontName', 'Arial', 'FontSize', options.AxisFontSize, 'LineWidth', 1.5);  % For tick labels
    hold on; grid on; axis equal;
    nIter = length(optimizationData.IterationData);

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

        % Gradient (global initial point) -> normalized
        gradientDirection = optimizationData.ProblemData(1).DesignSpaceLowerBound - ...
                            optimizationData.IterationData(1).OptimumCandidate;
        gradientDescent = -gradientDirection / (norm(gradientDirection) + eps);

        % Define arrow start and end
        gradientStartPoint = optimizationData.ProblemData(1).DesignSpaceLowerBound + 0.1;
        gradientEndPoint = gradientStartPoint + 0.25 * gradientDescent;

        hold on;
        h_quiver = quiver(gradientStartPoint(1), gradientStartPoint(2), ...
                          gradientEndPoint(1) - gradientStartPoint(1), ...
                          gradientEndPoint(2) - gradientStartPoint(2), ...
                          0, 'k', 'LineWidth', 2, 'MaxHeadSize', 2, ...
                          'DisplayName', 'Global Gradient');

        % Add legend
        legend(h_quiver, 'Location', 'northwest');
    end

    % Set axis limits
    xLimits = [min(optimizationData.ProblemData(1).DesignSpaceLowerBound(1)), ...
               max(optimizationData.ProblemData(1).DesignSpaceUpperBound(1))];
    yLimits = [min(optimizationData.ProblemData(1).DesignSpaceLowerBound(2)), ...
               max(optimizationData.ProblemData(1).DesignSpaceUpperBound(2))];
    axis([xLimits, yLimits]);

    % Labels and title
    titleHandle = title(sprintf('%s optimization progress - Iteration Index (k) = 0', options.OptimizationMethod), ...
                        'FontName', 'Arial', 'FontSize', options.TitleFontSize);
    xlabel(options.DesignVariableNames{1}, 'interpreter', 'latex', 'FontSize', options.AxisFontSize);
    ylabel(options.DesignVariableNames{2}, 'interpreter', 'latex', 'FontSize', options.AxisFontSize);

    % Initial point
    initialPoint = optimizationData.InitialData.OptimumCandidate;
    hInitial = plot(initialPoint(1), initialPoint(2), 'go', 'MarkerSize', 10, ...
                    'MarkerFaceColor', 'g', 'DisplayName', 'Initial Point');

    % Plot progress points and lines (all initially invisible)
    graphicsHandles3.points = gobjects(nIter, 1);
    graphicsHandles3.lines = gobjects(nIter-1, 1);
    currentPoint = initialPoint;

    for i = 1:nIter
        thisPoint = optimizationData.IterationData(i).OptimumCandidate;
        graphicsHandles3.points(i) = plot(thisPoint(1), thisPoint(2), 'ko', ...
            'MarkerSize', 10, 'MarkerFaceColor', 'k', 'Visible', 'off', ...
            'DisplayName', 'Intermediate point');

        if i > 1
            prevPoint = optimizationData.IterationData(i-1).OptimumCandidate;
            graphicsHandles3.lines(i-1) = plot([prevPoint(1), thisPoint(1)], ...
                                               [prevPoint(2), thisPoint(2)], ...
                                               'k--', 'LineWidth', 2, 'Visible', 'off', ...
                                               'DisplayName', 'Path');
        end
    end

    % Final point
    finalPoint = optimizationData.IterationData(nIter).OptimumCandidate;
    graphicsHandles3.points(nIter) = plot(finalPoint(1), finalPoint(2), 'o', ...
        'MarkerSize', 10, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', ...
        'DisplayName', 'Final point');

    % Slider setup
    sliderHandle = uicontrol('Style', 'slider', ...
        'Min', 0, 'Max', nIter, 'Value', 0, ...
        'SliderStep', [1/nIter, 1/nIter], ...
        'Units', 'normalized', ...
        'Position', [0.25, 0.01, 0.5, 0.04], ...
        'Callback', @sliderCallback);

    % Add dummy line to maintain consistent legend entries
    dummyPath = plot(nan, nan, 'k--', 'LineWidth', 2, 'DisplayName', 'Path');
    legend([hInitial, graphicsHandles3.points(nIter), dummyPath], ...
           'Location', 'northeast');

    % Save initial GIF frame
    frame = getframe(gcf);
    img = frame2im(frame);
    [A, map] = rgb2ind(img, 256);
    imwrite(A, map, gif_filename, 'gif', 'LoopCount', inf, 'DelayTime', frameDelay, 'WriteMode', 'overwrite');

    % Animate for gif
    for i = 1:nIter
        set(sliderHandle, 'Value', i);
        sliderCallback(sliderHandle);
        drawnow;

        frame = getframe(gcf);
        img = frame2im(frame);
        [A, map] = rgb2ind(img, 256);
        imwrite(A, map, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frameDelay);
    end

    % === Nested Callback Function ===
    function sliderCallback(src, ~)
        val = round(get(src, 'Value'));
        set(src, 'Value', val);
        titleHandle.String = sprintf('%s optimization progress - Iteration Index (k) = %d', options.OptimizationMethod, val);

        for j = 1:nIter
            set(graphicsHandles3.points(j), 'Visible', j <= val);
        end
        for j = 1:nIter-1
            set(graphicsHandles3.lines(j), 'Visible', j <= val-1);
        end
    end
end
