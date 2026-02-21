function plot_contour_objective_constraint_results(optimizationData, constraints, varargin)
%PLOT_CONTOUR_OBJECTIVE_CONSTRAINT_RESULTS Visualizes constrained optimization with animated progress.
%
%   PLOT_CONTOUR_OBJECTIVE_CONSTRAINT_RESULTS(OPTIMIZATIONDATA, CONSTRAINTS, ...)
%   creates an animated visualization of constrained optimization progress with:
%       - Contour plot of objective function (2D only)
%       - Constraint boundary visualization
%       - Sequential iteration path animation
%       - Automatic GIF/video export
%       - Real-time progress display
%
%   INPUTS:
%       OPTIMIZATIONDATA - Structure containing optimization history:
%           .ProblemData - Problem definition with bounds and objective function
%           .InitialData - Initial point information
%           .IterationData - Array of iteration records with OptimumCandidate
%       CONSTRAINTS - Function handle for constraint evaluation (returns constraint violations)
%
%   NAME-VALUE PAIR ARGUMENTS:
%       'SaveFolder' - Directory to save outputs (default: [])
%       'SaveFigureOptions' - Options for figure saving (default: {})
%       'CreateVideoGifIterations' - Enable animation export (default: true)
%       'FramesPerSecond' - Animation frame rate (default: 5)
%       'ImwriteOptions' - Options for GIF creation (default: {})
%       'VideoWriterProfile' - Video codec profile (default: 'Motion JPEG 2000')
%       'VideoWriterOptions' - Video writer settings (default: {'LosslessCompression', true})
%       'TitleFontSize' - Font size for title (default: 14)
%       'AxisFontSize' - Font size for axes (default: 12)
%       'DesignVariableNames' - Custom names for variables (default: x₁,x₂,...)
%       'UseDesignSpaceAsAxisLimits' - Use design space bounds for limits (default: true)
%       'PlotContour2DUseDesignSpaceLimit' - Use design space for contour (default: true)
%       'PlotContour2DOptions' - Additional contour plot options (default: {})
%       'ArrowScale' - Scaling factor for direction arrows (default: 0.01)
%
%   FUNCTIONALITY:
%       1. Generates contour plot of objective function (2D problems)
%       2. Visualizes constraint boundaries as red contour lines
%       3. Animates optimization path with sequential point plotting
%       4. Displays initial point (black), intermediate points (black), final point (blue)
%       5. Connects iterations with magenta dashed lines
%       6. Exports animation as both GIF and video files
%       7. Maintains design space aspect ratio and bounds
%
%   ANIMATION SEQUENCE:
%       - Plots objective function contour background
%       - Overlays constraint boundary (red line at constraint = 0)
%       - Shows initial point
%       - Sequentially reveals each iteration point
%       - Draws connection lines between consecutive points
%       - Highlights final optimum point
%       - Captures each frame for GIF/video export
%
%   CONSTRAINT VISUALIZATION:
%       - Evaluates constraint function over design space grid
%       - Plots constraint boundary as red contour line where constraint = 0
%       - Feasible region is where constraint ≤ 0
%       - Infeasible region is where constraint > 0
%
%   NOTES:
%       - Designed specifically for 2D constrained optimization problems
%       - Uses Arial font consistently throughout
%       - Automatically handles design space bounds
%       - Creates high-quality output for presentations/publications
%       - Animation timing controlled by FramesPerSecond parameter
%       - Particularly useful for penalty/barrier method visualization
%       - Constraint function must return scalar violation values
%
%   EXAMPLE:
%       % Visualize constrained optimization progress
%       constraintFunc = @(x) x(:,1) - 0.5;  % x1 >= 0.5 constraint
%       plot_contour_objective_constraint_results(optData, constraintFunc, ...
%           'SaveFolder', 'results/', 'FramesPerSecond', 8, ...
%           'CreateVideoGifIterations', true)
%
%   SEE ALSO:
%       plot_contour_objective_constraint_results_slider, plot_point_optimization_contour_2d,
%       optimization_barrier_penalty, optimization_sequential_quadratic_programing
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
    parser.addParameter('DesignVariableNames',{});
    parser.addParameter('UseDesignSpaceAsAxisLimits',true);
    parser.addParameter('PlotContour2DUseDesignSpaceLimit',true);
    parser.addParameter('PlotContour2DOptions',{});
    parser.addParameter('ArrowScale',0.01, @(x) isnumeric(x) && isscalar(x));
    parser.parse(varargin{:});
    options = parser.Results;

    % design variable names
    if(isempty(options.DesignVariableNames))
        for i=1:length(optimizationData.InitialData.OptimumCandidate)
            options.DesignVariableNames{i} = sprintf('$$x_{%d}$$',i);
        end
    end

    % Setup for video/gif saving
    defaultVideoWriterOptions = {'FrameRate', options.FramesPerSecond};
    [~, videoWriterOptions] = merge_name_value_pair_argument(defaultVideoWriterOptions, options.VideoWriterOptions);
    defaultImwriteOptions = {};
    [~, imwriteOptions] = merge_name_value_pair_argument(defaultImwriteOptions, options.ImwriteOptions);
    frameDelay = 1 / options.FramesPerSecond;

    %% Figure 1

    gif_filename = fullfile(options.SaveFolder, 'Optimization_progress.gif');
    
    % Create figure
    figureHandle = figure;
    set(gcf, 'Color', 'w');  % Set the figure background to white
    hold on;
    grid on;
    axis equal;

    nIter = length(optimizationData.IterationData);

    % Define grid for plotting
    [x1, x2] = meshgrid(linspace(optimizationData.ProblemData.DesignSpaceLowerBound(1), optimizationData.ProblemData.DesignSpaceUpperBound(1), 100),...
        linspace(optimizationData.ProblemData.DesignSpaceLowerBound(2), optimizationData.ProblemData.DesignSpaceUpperBound(2), 100));
    [Z,~] = evaluate_optimization_objective(optimizationData.ProblemData.ObjectiveFunction,[x1(:),x2(:)]);
    Z = reshape(Z, size(x1)); 

    % Plot contour once (for 2D only)
    if length(optimizationData.InitialData.OptimumCandidate) == 2
        if options.PlotContour2DUseDesignSpaceLimit
            contourDesignSpace = [optimizationData.ProblemData(1).DesignSpaceLowerBound;
                                  optimizationData.ProblemData(1).DesignSpaceUpperBound];
        else
            contourDesignSpace = [min(optimumCandidate,[],1); max(optimumCandidate,[],1)];
        end
        hContour = plot_point_optimization_contour_2d(figureHandle, optimizationData.ProblemData(1).ObjectiveFunction, contourDesignSpace, ...
            'DesignVariableNames', options.DesignVariableNames, ...
            'AxisFontSize', options.AxisFontSize, ...
            options.PlotContour2DOptions{:});

        % Hide all graphics objects created in this function from legend
        if isgraphics(hContour)
            set(findall(hContour, '-property', 'HandleVisibility'), 'HandleVisibility', 'off');
        end
    end

    % Plot Constraint contour (x1 - 0.5 ≥ 0 → x1 ≥ 0.5)
    % Evaluate constraint function over the grid
    X_flat = [x1(:), x2(:)];  % Convert grid to list of points
    [C, ~] = constraints(X_flat);  % Compute constraint violation
    C = reshape(C, size(x1));  % Reshape to match grid size

    hold on;

    
    % **1. Plot Constraint Boundary Using Contour**
    contour(x1, x2, C, [0 0], 'LineColor', 'r', 'LineWidth', 2, 'DisplayName', 'Constraint Boundary');

    % Extract problem bounds
    % Set axis limits
    xLower = optimizationData.ProblemData(1).DesignSpaceLowerBound(1);
    xUpper = optimizationData.ProblemData(1).DesignSpaceUpperBound(1);
    xLimits = [min(xLower), max(xUpper)];
    
    yLower = optimizationData.ProblemData(1).DesignSpaceLowerBound(2);
    yUpper = optimizationData.ProblemData(1).DesignSpaceUpperBound(2);
    yLimits = [min(yLower), max(yUpper)];
    if length(optimizationData.InitialData.OptimumCandidate) == 3
        zLimits = [optimizationData.ProblemData(3).DesignSpaceLowerBound, optimizationData.ProblemData(3).DesignSpaceUpperBound];
        axis([xLimits, yLimits, zLimits]);
    else
        axis([xLimits, yLimits]);
    end
    
    % Plot Format
    title('Penalty barrier optimization progess', ...
        'FontName', 'Arial', 'FontSize', options.TitleFontSize);
    xlabel(options.DesignVariableNames{1},'interpreter','latex','FontSize',options.AxisFontSize);
    ylabel(options.DesignVariableNames{2},'interpreter','latex','FontSize',options.AxisFontSize);
    if length(optimizationData.InitialData.OptimumCandidate) == 3
        zlabel(options.DesignVariableNames{3},'interpreter','latex','FontSize',options.AxisFontSize);
    end
    
    % Store previous points for line connections
    previousX = [];
    previousY = [];
    
    % Start video/gif
    if (options.CreateVideoGifIterations)
        filename = [options.SaveFolder, 'OptimizationProgression'];
        videoHandle = VideoWriter(filename, options.VideoWriterProfile);
        for i = 1:2:length(videoWriterOptions)
            videoHandle.(videoWriterOptions{i}) = videoWriterOptions{i+1};
        end
        open(videoHandle);
    end
    
    % Capture the first frame for GIF (after contour plot is created)
    frame = getframe(gcf);
    img = frame2im(frame);
    [A, map] = rgb2ind(img, 256);
    
    % Save as GIF (first frame)
    imwrite(A, map, filename, 'gif', 'LoopCount', inf, 'DelayTime', frameDelay, 'WriteMode', 'overwrite');
    
    % Pause before next iteration
    pause(frameDelay);

    % Plot initial point
    initialPoint = optimizationData.InitialData.OptimumCandidate;
    % plot(initialPoint(1), initialPoint(2), 'MarkerSize', 10, 'MarkerFaceColor', [0.5, 0.5, 0.5]);

    % Loop through iterations, displaying points one by one
    for i = 1:nIter-1        
        % First, plot the current point
        if i == 1
            % Extract current and next points
            currentPoint = initialPoint;

            plot(initialPoint(1), initialPoint(2), 'o', ... % specify marker type
            'MarkerSize', 10, ...
            'MarkerEdgeColor', 'k', ... % black edge for visibility
            'MarkerFaceColor', 'k', ...
            'DisplayName', 'initial point');

        else
            % Extract current and next points
            currentPoint = optimizationData.IterationData(i).OptimumCandidate;

            plot(currentPoint(1), currentPoint(2), 'ko', ...
                 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', 'Current point'); % black for current point
        end

        % Store previous points for connection
        previousX = [previousX, currentPoint(1)];
        previousY = [previousY, currentPoint(2)];
        
        % Pause to make the appearance gradual
        pause(frameDelay);

        nextPoint = optimizationData.IterationData(i+1).OptimumCandidate;
        
        % After the pause, plot the connection line (magenta) between the current and next point
        plot([previousX(end), nextPoint(1)], [previousY(end), nextPoint(2)], 'm--', 'LineWidth', 2);

        % Update axes
        set(gca, 'FontSize', 12, 'LineWidth', 1.5);
        xlim(xLimits); ylim(yLimits);
        if length(optimizationData.InitialData.OptimumCandidate) == 3
            zlim(zLimits);
        end
    
        drawnow;
        pause(frameDelay);

        % Save video/GIF frame
        if options.CreateVideoGifIterations
            currentFrame = getframe(gcf);
            writeVideo(videoHandle, currentFrame);
    
            [A, map] = rgb2ind(frame2im(currentFrame), 256);
            if i == 1
                imwrite(A, map, [filename, '.gif'], 'gif', 'LoopCount', inf, 'DelayTime', frameDelay, imwriteOptions{:});
            else
                imwrite(A, map, [filename, '.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', frameDelay, imwriteOptions{:});
            end
        end

    end
    
    % Plot the final optimum point in blue
    finalPoint = optimizationData.IterationData(nIter).OptimumCandidate;
    plot(finalPoint(1), finalPoint(2), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b'); % Blue for final optimum
        
    % Capture the final frame for GIF
    if options.CreateVideoGifIterations
        currentFrame = getframe(gcf);
        writeVideo(videoHandle, currentFrame);

        [A, map] = rgb2ind(frame2im(currentFrame), 256);
        if i == 1
            imwrite(A, map, [filename, '.gif'], 'gif', 'LoopCount', inf, 'DelayTime', frameDelay, imwriteOptions{:});
        else
            imwrite(A, map, [filename, '.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', frameDelay, imwriteOptions{:});
        end
    end
    
    % Close video handle
    if (options.CreateVideoGifIterations)
        close(videoHandle);
    end
end