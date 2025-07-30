function plot_point_gradient_progress_slider(optimizationData, objectiveFunction, varargin)
%PLOT_POINT_GRADIENT_PROGRESS_SLIDER Visualize gradient-based optimization progress with interactive controls.
%
%   PLOT_POINT_GRADIENT_PROGRESS_SLIDER(OPTIMIZATIONDATA, OBJECTIVEFUNCTION, ...)
%   creates an interactive visualization of gradient-based optimization progress,
%   including:
%       - Main optimization path with gradient/search direction arrows
%       - Line search process animation
%       - Final optimization progress overview
%
%   INPUTS:
%       OPTIMIZATIONDATA     - Structure containing optimization history:
%           .ProblemData     - Problem definition with bounds
%           .InitialData     - Initial point information
%           .IterationData   - Array of iteration records
%       OBJECTIVEFUNCTION    - Function handle for objective evaluation
%
%   NAME-VALUE PAIR ARGUMENTS:
%       'OptimizationMethod' - Name of method used (default: 'Conjugate gradient')
%       'SaveFolder'         - Directory to save outputs (default: [])
%       'CreateVideoGifIterations' - Enable video/GIF creation (default: true)
%       'FramesPerSecond'    - Animation frame rate (default: 5)
%       'VideoWriterProfile' - Video encoding profile (default: 'Motion JPEG 2000')
%       'DesignVariableNames'- Custom names for variables (default: x₁,x₂,...)
%       'ArrowScale'         - Scaling factor for direction arrows (default: 0.01)
%       (Complete list of 15+ parameters available in function)
%
%   OUTPUTS:
%       None (generates figures and optional video/GIF files)
%
%   FUNCTIONALITY:
%       1. Creates three integrated visualizations:
%          a) Main optimization path with interactive slider
%          b) Line search process animation per iteration
%          c) Final optimization progress overview
%       2. Supports both 2D and 3D problems
%       3. Generates:
%          - Contour plots (2D)
%          - Gradient/search direction vectors
%          - Iteration path with markers
%          - Line search bracketing visualization
%       4. Optional video/GIF output of entire process
%
%   NOTES:
%       - Uses Arial font consistently throughout
%       - Maintains design space aspect ratio
%       - Automatically handles both bracketing and sectioning line searches
%       - Slider allows interactive exploration of iterations
%
%   DEPENDENCIES:
%       - plot_point_optimization_contour_2d.m
%       - evaluate_optimization_objective.m
%       - merge_name_value_pair_argument.m
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor, Contributor)
%   Copyright 2025 Ali Abbas Kapadia (Main Author)
%   SPDX-License-Identifier: Apache-2.0

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
%    parser.addParameter('VideoWriterOptions', {});
    parser.parse(varargin{:});
    options = parser.Results;

    set(groot, ...
    'DefaultAxesFontName', 'Arial', ...
    'DefaultAxesFontSize', options.AxisFontSize);

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
    
    % Create figure
    figureHandle = figure;
    set(gcf, 'Color', 'w');
    hold on;
    axis equal;
    grid on;

    % Set axis limits
    xLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(1), optimizationData.ProblemData.DesignSpaceUpperBound(1)];
    yLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(2), optimizationData.ProblemData.DesignSpaceUpperBound(2)];
    if length(optimizationData.InitialData.OptimumCandidate) == 3
        zLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(3), optimizationData.ProblemData.DesignSpaceUpperBound(3)];
        axis([xLimits, yLimits, zLimits]);
    else
        axis([xLimits, yLimits]);
    end

    % Plot contour once (2D only)
    if length(optimizationData.InitialData.OptimumCandidate) == 2
        if options.PlotContour2DUseDesignSpaceLimit
            contourDesignSpace = [optimizationData.ProblemData(1).DesignSpaceLowerBound;
                                  optimizationData.ProblemData(1).DesignSpaceUpperBound];
        else
            contourDesignSpace = [min(optimumCandidate, [], 1); max(optimumCandidate, [], 1)];
        end

        hContour = plot_point_optimization_contour_2d(figureHandle, ...
            optimizationData.ProblemData(1).ObjectiveFunction, ...
            contourDesignSpace, ...
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

    xlabel(options.DesignVariableNames{1}, 'interpreter','latex','FontSize', options.AxisFontSize);
    ylabel(options.DesignVariableNames{2}, 'interpreter','latex','FontSize', options.AxisFontSize);
    if length(optimizationData.InitialData.OptimumCandidate) == 3
        zlabel(options.DesignVariableNames{3}, 'interpreter','latex','FontSize', options.AxisFontSize);
    end

    % Video settings
    if options.CreateVideoGifIterations
        filename = [options.SaveFolder, 'IterationProgression'];
        videoHandle = VideoWriter(filename, options.VideoWriterProfile);
        for i = 1:2:length(videoWriterOptions)
            videoHandle.(videoWriterOptions{i}) = videoWriterOptions{i+1};
        end
        open(videoHandle);
    end

    % Iteration setup
    nIter = length(optimizationData.IterationData);
    graphicsHandles1 = struct('pointPrev', [], 'pointNext', [], 'arrowGrad', [], 'arrowSearch', [], 'stepLine', []);
    fixedArrowSize = options.ArrowScale;

    % Plot all iterations first (for video/GIF)
    for i = 1:nIter
        update_plot(i);
        title(sprintf('%s optimization process - Iteration %d', options.OptimizationMethod, i), ...
            'FontName', 'Arial', 'FontSize', options.TitleFontSize);
        drawnow;
        pause(frameDelay);

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

    if options.CreateVideoGifIterations
        close(videoHandle);
    end

    % Add UI slider for interactive iteration view
    if nIter > 1
        sliderHandle = uicontrol('Style', 'slider', 'Min', 1, 'Max', nIter, 'Value', nIter, ...
            'SliderStep', [1/(nIter-1) , 1/(nIter-1)], ...
            'Units', 'normalized', 'Position', [0.25 0.01 0.5 0.03]);

        addlistener(sliderHandle, 'Value', 'PostSet', @(src, event) ...
            update_plot(round(get(sliderHandle, 'Value'))));
    end

    % Function to update plots at a given iteration index
    function update_plot(i)
        isFinalPoint = (i == nIter);
        next = optimizationData.IterationData(i).OptimumCandidate;
        grad = optimizationData.IterationData(i).Gradient;
        search = optimizationData.IterationData(i).SearchDirection;

        if(i==1)
            current = optimizationData.InitialData.OptimumCandidate;
        else
            current = optimizationData.IterationData(i-1).OptimumCandidate;
        end

        % Delete old elements
        names = fieldnames(graphicsHandles1);
        for k = 1:length(names)
            h = graphicsHandles1.(names{k});
            if ~isempty(h)
                validHandles = h(isgraphics(h));
                delete(validHandles);
            end
        end

        if length(current) == 2
            % Current point (black)
            graphicsHandles1.pointPrev = plot(current(1), current(2), 'o', ...
                                 'MarkerSize', 10, 'MarkerFaceColor', [0.5 0.5 0.5],'MarkerEdgeColor', [0.5 0.5 0.5], 'DisplayName', 'Previous point');

            % Gradient arrow (red)
            graphicsHandles1.arrowGrad = quiver(current(1), current(2), grad(1)*fixedArrowSize, grad(2)*fixedArrowSize, ...
                                   0, 'Color', 'r', 'LineWidth', 1.5, 'MaxHeadSize', 2, 'DisplayName', 'Gradient');

            % Search direction (dark blue)
            graphicsHandles1.arrowSearch = quiver(current(1), current(2), search(1)*fixedArrowSize, search(2)*fixedArrowSize, ...
                                     0, 'Color', [0 0 0.5], 'LineWidth', 1.5, 'MaxHeadSize', 2, 'DisplayName', 'Search direction');

            % Next point (grey)
            graphicsHandles1.pointNext = plot(next(1), next(2), 'ko', ...
                                 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', ...
                                 'DisplayName', 'Current point');

            % Line between points
            graphicsHandles1.stepLine = plot([current(1), next(1)], [current(2), next(2)], 'k--', ...
                                'LineWidth', 2, 'HandleVisibility', 'off');

            % If this is the final point, plot in yellow with black border
            if isFinalPoint
                graphicsHandles1.finalPoint = plot(next(1), next(2), 'o', ...
                                     'MarkerSize', 10, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', ...
                                     'DisplayName', 'Final point');
            end

            xlim(xLimits); ylim(yLimits);
        else
            graphicsHandles1.pointNext = plot3(current(1), current(2), current(3), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', 'Previous point');
            graphicsHandles1.arrowGrad = quiver3(current(1), current(2), current(3), grad(1)*fixedArrowSize, grad(2)*fixedArrowSize, grad(3)*fixedArrowSize, ...
                                                0, 'Color', 'r', 'LineWidth', 1.5, 'MaxHeadSize', 2, 'DisplayName', 'Gradient');
            graphicsHandles1.arrowSearch = quiver3(current(1), current(2), current(3), search(1)*fixedArrowSize, search(2)*fixedArrowSize, search(3)*fixedArrowSize, ...
                                                  0, 'Color', [0 0 0.5], 'LineWidth', 1.5, 'MaxHeadSize', 2, 'DisplayName', 'Search direction');
            graphicsHandles1.pointPrev = plot3(next(1), next(2), next(3), 'o', ...
                                 'MarkerSize', 10, 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], ...
                                 'DisplayName', 'Current point');
            graphicsHandles1.stepLine = plot3([current(1), next(1)], [current(2), next(2)], [current(3), next(3)], 'k--', 'LineWidth', 2);

            % If this is the final point, plot in yellow with black border
            if isFinalPoint
                graphicsHandles1.finalPoint = plot(next(1), next(2), next(3), 'o', ...
                                     'MarkerSize', 10, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', ...
                                     'DisplayName', 'Final point');
            end

            xlim(xLimits); ylim(yLimits); zlim(zLimits);
        end

        title(sprintf('%s optimization process - k = %d', options.OptimizationMethod, i), ...
            'FontName', 'Arial', 'FontSize', options.TitleFontSize);
        set(gca, 'FontName', 'Arial', 'FontSize', options.AxisFontSize, 'LineWidth', 1.5);
        legend('Location', 'northeast', 'FontSize', 12);
    end
    %% Figure 2: Line Search Animation

    % Create a LineSearch subfolder
    lineSearchFolder = fullfile(options.SaveFolder, 'LineSearch');
    if ~exist(lineSearchFolder, 'dir')
        mkdir(lineSearchFolder);
    end

    % Create figure
    figure;
    set(gcf, 'Color', 'w');  % Set the figure background to white
    axis equal;
    grid on;

    % Start video/gif
    if (options.CreateVideoGifIterations)
        filename = fullfile(lineSearchFolder, 'StepSizeProgression');
        videoHandle = VideoWriter(filename, options.VideoWriterProfile);
        for i = 1:2:length(videoWriterOptions)
            videoHandle.(videoWriterOptions{i}) = videoWriterOptions{i+1};
        end
        open(videoHandle);
    end

    if isequal(optimizationData.ProblemData.Options.LineSearchFunction, @line_search_backtracking)

        % Iterate over each optimization step
        for i = 1:nIter
            % Extract current design and gradient
            if(i==1)    
                currentDesign   = optimizationData.InitialData.OptimumCandidate;
            else
                currentDesign   = optimizationData.IterationData(i-1).OptimumCandidate;
            end
            searchDirection = optimizationData.IterationData(i).SearchDirection;
            objectiveValue  = optimizationData.IterationData(i).StepSizeData.ObjectiveValues;

            % Get step sizes used in the iteration
            stepsize_history = optimizationData.IterationData(i).StepSizeData.StepSizes; % Extracted step sizes
            final_stepsize = stepsize_history(end); % Last step size is α*

            % Define step-size sweep
            stepSize_sweep = linspace(0, max(stepsize_history)+0.2, 100); % Sweeping from 0 to 1
            nextDesign_sweep = currentDesign + stepSize_sweep' * searchDirection;

            [objective_sweep,~] = evaluate_optimization_objective(objectiveFunction, nextDesign_sweep);

            % Plot the line search (Objective function curve)
            clf; hold on;
            plot(stepSize_sweep, objective_sweep, 'b-', 'LineWidth', 2); % Plot objective curve
            xlabel('Step size ($\alpha$)', 'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
            ylabel('Objective Value ($f(\mathbf{x})$)',  'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
            title(['k = ', num2str(i)], 'FontName', 'Arial', 'FontSize', options.TitleFontSize);
            set(gca, 'FontName', 'Arial', 'FontSize', options.AxisFontSize, 'LineWidth', 1.5);  % For tick labels
            grid on;

            % Initialize filename for this iteration's GIF
            gif_filename = fullfile(lineSearchFolder, sprintf('LineSearch_Iteration_%d.gif', i));

            % Plot step-size intersections with different colors
            for j = 1:length(stepsize_history)

                % Determine marker style
                if j == 1
                    markerColor = [0.5, 0.5, 0.5]; % Gray for first intersection
                    markerSize = 8;
                elseif j == length(stepsize_history)
                    markerColor = 'b'; % Blue for final intersection
                    markerSize = 10;
                else
                    markerColor = 'r'; % Red for intermediate intersections
                    markerSize = 6;
                end

                % Plot the intersection points
                plot(stepsize_history(j), objectiveValue(j), 'o', ...
                     'MarkerSize', markerSize, 'MarkerEdgeColor', 'k', ...
                     'MarkerFaceColor', markerColor, 'LineWidth', 1.5);

                % Plot the vertical line for this step-size intersection
                xline(stepsize_history(j), '--k', ['\alpha_', num2str(j)], ...
                      'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'right', 'FontSize', 10);

                hold off;

                pause(frameDelay); % Pause for GIF effect

                % Save video/GIF frame
                if (options.CreateVideoGifIterations)
                    currentFrame = getframe(gcf);
                    writeVideo(videoHandle, currentFrame);

                    [A, map] = rgb2ind(frame2im(getframe(gcf)), 256); % Convert to indexed image
                    if (j == 1)
                        % For the first frame, create the GIF for this iteration
                        imwrite(A, map, gif_filename, 'gif', 'LoopCount', inf, 'DelayTime', frameDelay);
                    else
                        % Append subsequent frames to the GIF for this iteration
                        imwrite(A, map, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frameDelay);
                    end
                end
            end

            % Plot final step size as a red dotted line
            xline(final_stepsize, '--r', '\alpha_*', 'LineWidth', 2, ...
                  'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left', 'FontSize', 12);

            % Save the final video/GIF frame for this iteration
            if (options.CreateVideoGifIterations)
                currentFrame = getframe(gcf);
                writeVideo(videoHandle, currentFrame);

                [A, map] = rgb2ind(frame2im(getframe(gcf)), 256);
                imwrite(A, map, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frameDelay);
            end

            % Pause before the next iteration
            pause(frameDelay);
        end

    elseif isequal(optimizationData.ProblemData.Options.LineSearchFunction, @line_search_golden_ratio)
        for i = 1:nIter
            % Extract current design and gradient
            if(i==1)
                currentDesign   = optimizationData.InitialData.OptimumCandidate;
            else
                currentDesign   = optimizationData.IterationData(i-1).OptimumCandidate;
            end
            searchDirection = optimizationData.IterationData(i).SearchDirection;
            objectiveValue  = optimizationData.IterationData(i).StepSizeData.ObjectiveValues;

            % Get step sizes used in the iteration
            stepsize_history = optimizationData.IterationData(i).StepSizeData.StepSizes; % Extracted step sizes
            final_stepsize = stepsize_history(end); % Last step size is α*

            % Define step-size sweep
            stepSize_sweep = linspace(0, max(stepsize_history)+0.2, 100); % Sweeping from 0 to 1
            nextDesign_sweep = currentDesign + stepSize_sweep' * searchDirection;

            [objective_sweep,~] = evaluate_optimization_objective(objectiveFunction, nextDesign_sweep);

            % Plot the line search (Objective function curve)
            clf; hold on;
            plot(stepSize_sweep, objective_sweep, 'b-', 'LineWidth', 2); % Plot objective curve
            xlabel('Step size ($\alpha$)', 'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
            ylabel('Objective Value ($f(\mathbf{x})$)',  'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
            title(['k = ', num2str(i)], 'FontName', 'Arial', 'FontSize', options.TitleFontSize);
            set(gca, 'FontName', 'Arial', 'FontSize', options.AxisFontSize, 'LineWidth', 1.5);  % For tick labels
            grid on;

            % Initialize filename for this iteration's GIF
            gif_filename = fullfile(lineSearchFolder, sprintf('LineSearch_Iteration_%d.gif', i));

            % Plot step-size intersections with different colors
            for j = 1:length(stepsize_history)

                % Sectioning in red leading to final point
                if j < length(stepsize_history)
                    if j <= 2
                        markerColor = [0.5, 0.5, 0.5]; % Grey for bracketing
                        markerSize = 6;
                    else
                        markerColor = 'r'; % Red for sectioning
                        markerSize = 6;
                    end
                else
                    markerColor = 'b'; % Final point in blue
                    markerSize = 10;
                end

                % Plot the intersection points
                plot(stepsize_history(j), objectiveValue(j), 'o', ...
                     'MarkerSize', markerSize, 'MarkerEdgeColor', 'k', ...
                     'MarkerFaceColor', markerColor, 'LineWidth', 1.5);

                % Plot the vertical line for this step-size intersection
                xline(stepsize_history(j), '--k', ['\alpha_', num2str(j)], ...
                      'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'right', 'FontSize', 10);

                % Save video/GIF frame after each point is plotted
                if options.CreateVideoGifIterations
                    currentFrame = getframe(gcf);
                    [A, map] = rgb2ind(frame2im(currentFrame), 256); % Convert to indexed image
                    if (j == 1)
                        % For the first frame, create the GIF for this iteration
                        imwrite(A, map, gif_filename, 'gif', 'LoopCount', inf, 'DelayTime', frameDelay);
                    else
                        % Append subsequent frames to the GIF for this iteration
                        imwrite(A, map, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frameDelay);
                    end
                end
                pause(frameDelay); % Pause for GIF effect (controls speed)
            end

            % Plot final step size as a red dotted line (or blue, based on your logic)
            xline(final_stepsize, '--r', '\alpha_*', 'LineWidth', 2, ...
                  'LabelVerticalAlignment', 'middle', 'LabelHorizontalAlignment', 'left', 'FontSize', 12);

            % Save the final video/GIF frame for this iteration
            if options.CreateVideoGifIterations
                currentFrame = getframe(gcf);
                [A, map] = rgb2ind(frame2im(currentFrame), 256);
                imwrite(A, map, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frameDelay);
            end

            % Pause before the next iteration
            pause(frameDelay);
        end
    end

    hold off

    %% Figure 3: Visualization of Optimization Progress
    
    gif_filename = fullfile(options.SaveFolder, 'Optimization_progress.gif');
    
    % Create figure
    figureHandleOptimization = figure;
    set(gcf, 'Color', 'w');  % Set the figure background to white\
    set(gca, 'FontName', 'Arial', 'FontSize', options.AxisFontSize, 'LineWidth', 1.5);
    hold on;
    grid on;
    axis equal;
    
    % Extract problem bounds
    xLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(1), optimizationData.ProblemData.DesignSpaceUpperBound(1)];
    yLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(2), optimizationData.ProblemData.DesignSpaceUpperBound(2)];
    
    % Plot contour once (2D only)
    if length(optimizationData.InitialData.OptimumCandidate) == 2
        if options.PlotContour2DUseDesignSpaceLimit
            contourDesignSpace = [optimizationData.ProblemData(1).DesignSpaceLowerBound;
                                  optimizationData.ProblemData(1).DesignSpaceUpperBound];
        else
            contourDesignSpace = [min(optimumCandidate, [], 1); max(optimumCandidate, [], 1)];
        end
    
        hContour = plot_point_optimization_contour_2d(figureHandleOptimization, ...
            optimizationData.ProblemData(1).ObjectiveFunction, ...
            contourDesignSpace, ...
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
    
    % Plot Format
    titleHandle = title(sprintf('%s optimization progress - k = 0', options.OptimizationMethod), ...
        'FontName', 'Arial', 'FontSize', options.TitleFontSize);
    xlabel(options.DesignVariableNames{1}, 'interpreter', 'latex', 'FontSize', options.AxisFontSize);
    ylabel(options.DesignVariableNames{2}, 'interpreter', 'latex', 'FontSize', options.AxisFontSize);
    if length(optimizationData.InitialData.OptimumCandidate) == 3
        zlabel(options.DesignVariableNames{3}, 'interpreter', 'latex', 'FontSize', options.AxisFontSize);
    end
    
    % Plot Initial Point (Green)
    initialPoint = optimizationData.InitialData.OptimumCandidate;
    hInitial = plot(initialPoint(1), initialPoint(2), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Initial Point');
    
    % Graphics handle containers
    graphicsHandles3.points = gobjects(nIter, 1);
    graphicsHandles3.lines = gobjects(nIter-1, 1);
    
    % Plot all points and lines (but make them invisible initially)
    for i = 1:nIter-1
        if(i==1)
            currentPoint = optimizationData.InitialData.OptimumCandidate;
        else
            currentPoint = optimizationData.IterationData(i-1).OptimumCandidate;
        end
        nextPoint = optimizationData.IterationData(i).OptimumCandidate;
    
        graphicsHandles3.points(i) = plot(currentPoint(1), currentPoint(2), 'ko', ...
            'MarkerSize', 10, 'MarkerFaceColor', 'k', 'Visible', 'off');
    
        graphicsHandles3.lines(i) = plot([currentPoint(1), nextPoint(1)], ...
                                        [currentPoint(2), nextPoint(2)], ...
                                        'k--', 'LineWidth', 2, 'Visible', 'off', 'DisplayName', 'Path');
    end
    
    % Final point
    finalPoint = optimizationData.IterationData(nIter).OptimumCandidate;
    graphicsHandles3.points(nIter) = plot(finalPoint(1), finalPoint(2), 'o', ...
        'MarkerSize', 10, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', ...
        'DisplayName', 'Final point');
    
    % Setup slider
    sliderHandle = uicontrol('Style', 'slider', ...
        'Min', 0, 'Max', nIter, 'Value', 0, ...
        'SliderStep', [1/(nIter), 1/(nIter)], ...
        'Units', 'normalized', ...
        'Position', [0.25, 0.01, 0.5, 0.04], ...
        'Callback', @sliderCallback);
    
    % Legend (create dummy for path if none drawn yet)
    dummyPath = plot(nan, nan, 'k--', 'LineWidth', 2, 'DisplayName', 'Path');
    legend([hInitial, graphicsHandles3.points(nIter), dummyPath], ...
           'Location', 'northeast');
    
    % Slider callback function
    function sliderCallback(src, ~)
        val = round(get(src, 'Value'));
        set(src, 'Value', val);
    
        % Update title
        titleHandle.String = sprintf('%s optimization progress - k = %d', options.OptimizationMethod, val);
    
        % Update visibility of points and lines
        for j = 1:nIter
            if j <= val
                set(graphicsHandles3.points(j), 'Visible', 'on');
            else
                set(graphicsHandles3.points(j), 'Visible', 'off');
            end
        end
        for j = 1:nIter-1
            if j <= val-1
                set(graphicsHandles3.lines(j), 'Visible', 'on');
            else
                set(graphicsHandles3.lines(j), 'Visible', 'off');
            end
        end
    end
    
    % Save initial frame of the GIF
    frame = getframe(gcf);
    img = frame2im(frame);
    [A, map] = rgb2ind(img, 256);
    imwrite(A, map, gif_filename, 'gif', 'LoopCount', inf, 'DelayTime', frameDelay, 'WriteMode', 'overwrite');
    
    % Optional: create gif frames for all iterations
    for i = 1:nIter
        set(sliderHandle, 'Value', i);
        sliderCallback(sliderHandle);
        drawnow;
        
        % Save frame
        frame = getframe(gcf);
        img = frame2im(frame);
        [A, map] = rgb2ind(img, 256);
        imwrite(A, map, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frameDelay);
    end

end
