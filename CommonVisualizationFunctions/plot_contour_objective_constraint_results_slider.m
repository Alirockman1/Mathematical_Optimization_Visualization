function plot_contour_objective_constraint_results_slider(optimizationData, constraints, varargin)
% PLOT_CONTOUR_OBJECTIVE_CONSTRAINT_RESULTS_SLIDER Visualizes 2D optimization with interactive slider
%
%   PLOT_CONTOUR_OBJECTIVE_CONSTRAINT_RESULTS_SLIDER(OPTIMIZATIONDATA, CONSTRAINTS)
%   creates an interactive visualization of 2D optimization progress with:
%     - Contour plot of objective function
%     - Constraint boundaries
%     - Optimization path animation
%     - Interactive slider to scrub through iterations
%     - Option to save as GIF/video
%
%   INPUTS:
%       optimizationData - Structure containing optimization results with fields:
%           .ProblemData - Problem definition including bounds and objective function
%           .IterationData - Array of structs with optimization history
%           .InitialData - Initial conditions and values
%       constraints - Function handle for constraint evaluation (returns constraint violations)
%
%   OPTIONAL PARAMETERS (name-value pairs):
%       'SaveFolder' - Path to save output figures/videos (default: [])
%       'SaveFigureOptions' - Options for figure saving (cell array)
%       'CreateVideoGifIterations' - Enable GIF/video creation (default: true)
%       'FramesPerSecond' - Animation frame rate (default: 5)
%       'ImwriteOptions' - Options for GIF creation (cell array)
%       'VideoWriterProfile' - Video profile for VideoWriter (default: 'Motion JPEG 2000')
%       'VideoWriterOptions' - Video writing options (cell array)
%       'TitleFontSize' - Font size for title (default: 14)
%       'AxisFontSize' - Font size for axis labels (default: 12)
%       'DesignVariableNames' - Custom names for design variables (cell array)
%       'UseDesignSpaceAsAxisLimits' - Use problem bounds for axes (default: true)
%       'PlotContour2DUseDesignSpaceLimit' - Use bounds for contour (default: true)
%       'PlotContour2DOptions' - Additional contour plot options (cell array)
%       'ArrowScale' - Scaling factor for gradient arrow (default: 0.01)
%
%   OUTPUTS:
%       Interactive figure showing:
%       1. Objective function contour plot
%       2. Constraint boundaries (red line)
%       3. Optimization path with:
%          - Initial point (green)
%          - Intermediate points (black)
%          - Final point (yellow)
%       4. Slider to navigate iterations
%       5. Optionally saves animation as GIF/video
%
%   KEY FEATURES:
%       - Interactive exploration of optimization history
%       - Clear visualization of objective/constraint landscape
%       - Professional formatting with LaTeX labels
%       - Comprehensive animation controls
%       - Flexible saving options
%
%   SEE ALSO:
%       plot_point_optimization_contour_2d, optimization_sequential_quadratic_programming
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

    % Parse input options
    parser = inputParser;
    parser.addParameter('SaveFolder', []);
    parser.addParameter('SaveFigureOptions', {});
    parser.addParameter('CreateVideoGifIterations', true);
    parser.addParameter('FramesPerSecond', 5);
    parser.addParameter('ImwriteOptions', {});
    parser.addParameter('VideoWriterProfile', 'Motion JPEG 2000');
    parser.addParameter('VideoWriterOptions', {'LosslessCompression', true});
    parser.addParameter('TitleFontSize', 14, @(x) isnumeric(x) && isscalar(x));
    parser.addParameter('AxisFontSize', 12, @(x) isnumeric(x) && isscalar(x));
    parser.addParameter('DesignVariableNames', {});
    parser.addParameter('UseDesignSpaceAsAxisLimits', true);
    parser.addParameter('PlotContour2DUseDesignSpaceLimit', true);
    parser.addParameter('PlotContour2DOptions', {});
    parser.addParameter('ArrowScale', 0.01, @(x) isnumeric(x) && isscalar(x));
    parser.parse(varargin{:});
    options = parser.Results;

    set(groot, ...
    'DefaultAxesFontName', 'Arial', ...
    'DefaultAxesFontSize', options.AxisFontSize);

    % Set up design variable names
    if isempty(options.DesignVariableNames)
        for i = 1:length(optimizationData.InitialData.OptimumCandidate)
            options.DesignVariableNames{i} = sprintf('$$x_{%d}$$', i);
        end
    end

    % Set up video/gif saving
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
        
        % Define grid for plotting
        [x1, x2] = meshgrid(linspace(optimizationData.ProblemData.DesignSpaceLowerBound(1), optimizationData.ProblemData.DesignSpaceUpperBound(1), 100),...
            linspace(optimizationData.ProblemData.DesignSpaceLowerBound(2), optimizationData.ProblemData.DesignSpaceUpperBound(2), 100));
        [Z,~] = evaluate_optimization_objective(optimizationData.ProblemData.ObjectiveFunction,[x1(:),x2(:)]);
        Z = reshape(Z, size(x1)); 
        
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
        hQuiver = quiver(gradientStartPoint(1), gradientStartPoint(2), ...
                          gradientEndPoint(1) - gradientStartPoint(1), ...
                          gradientEndPoint(2) - gradientStartPoint(2), ...
                          0, 'k', 'LineWidth', 2, 'MaxHeadSize', 2, ...
                          'DisplayName', 'Global Gradient');

        % Add legend
        legend(hQuiver, 'Location', 'northwest');
    end

    % Plot Constraint contour (x1 - 0.5 ≥ 0 → x1 ≥ 0.5)
    % Evaluate constraint function over the grid
    xGrid = [x1(:), x2(:)];                                                % Convert grid to list of points
    [C, ~] = constraints(xGrid);                                           % Compute constraint violation
    C = reshape(C, size(x1));                                              % Reshape to match grid size

    hold on;

    % **1. Plot Constraint Boundary Using Contour**
    contour(x1, x2, C, [0 0], 'LineColor', 'r', 'LineWidth', 2, 'DisplayName', 'Constraint Boundary');

    % Set axis limits
    xLimits = [min(optimizationData.ProblemData(1).DesignSpaceLowerBound(1)), ...
               max(optimizationData.ProblemData(1).DesignSpaceUpperBound(1))];
    yLimits = [min(optimizationData.ProblemData(1).DesignSpaceLowerBound(2)), ...
               max(optimizationData.ProblemData(1).DesignSpaceUpperBound(2))];
    axis([xLimits, yLimits]);

    % Labels and title
    titleHandle = title('Penalty barrier optimization - Iteration Index (k) = 0', ...
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
        nextPoint = optimizationData.IterationData(i).OptimumCandidate;

        graphicsHandles3.points(i) = plot(nextPoint(1), nextPoint(2), 'ko', ...
            'MarkerSize', 10, 'MarkerFaceColor', 'k', 'Visible', 'off', ...
            'DisplayName', 'Intermediate point');

        if i == 1
            % Create line from initial to first iteration
            graphicsHandles3.lines(1) = plot([initialPoint(1), nextPoint(1)], ...
                                             [initialPoint(2), nextPoint(2)], ...
                                             'k--', 'LineWidth', 2, 'Visible', 'off', ...
                                             'DisplayName', 'Path');
        elseif i > 1
            % Create line from previous iteration to current
            graphicsHandles3.lines(i) = plot([currentPoint(1), nextPoint(1)], ...
                                             [currentPoint(2), nextPoint(2)], ...
                                             'k--', 'LineWidth', 2, 'Visible', 'off', ...
                                             'DisplayName', 'Path');
        end

        currentPoint = nextPoint;
    end

    % Final point (only visible at last iteration)
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
        titleHandle.String = sprintf('Penalty barrier optimization - Iteration Index (k) = %d', val);

        % Update visibility of points and lines based on the current slider value
        for j = 1:nIter
            set(graphicsHandles3.points(j), 'Visible', j <= val); 
        end
        for j = 1:val
            if isgraphics(graphicsHandles3.lines(j))
                set(graphicsHandles3.lines(j), 'Visible', 'on');
            end
        end
        for j = val+1:length(graphicsHandles3.lines)
            if isgraphics(graphicsHandles3.lines(j))
                set(graphicsHandles3.lines(j), 'Visible', 'off');
            end
        end

        % Final point visibility (only show at last iteration)
        if val == nIter
            set(graphicsHandles3.points(nIter), 'Visible', 'on');
        else
            set(graphicsHandles3.points(nIter), 'Visible', 'off');
        end
    end
end

