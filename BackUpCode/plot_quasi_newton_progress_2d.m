function plot_quasi_newton_progress_2d(optimizationData, objectiveFunction, varargin)
    parser = inputParser;
    parser.addParameter('SaveFolder', []);
    parser.addParameter('SaveFigureOptions', {});
    parser.addParameter('CreateVideoGifIterations', true);
    parser.addParameter('FramesPerSecond', 5);
    parser.addParameter('ImwriteOptions', {});
    parser.addParameter('VideoWriterProfile', 'Motion JPEG 2000');
    parser.addParameter('VideoWriterOptions', {'LosslessCompression', true});
    parser.addParameter('LineSearchFunction', @line_search_golden_ratio, @(x) isa(x, 'function_handle'));
    parser.addParameter('AxisFontSize', 12, @(x) isnumeric(x) && isscalar(x));
    parser.addParameter('DesignVariableNames',{});
    parser.addParameter('UseDesignSpaceAsAxisLimits',true);
    parser.addParameter('PlotContour2DUseDesignSpaceLimit',true);
    parser.addParameter('PlotContour2DOptions',{});
    parser.addParameter('ArrowScale',0.01, @(x) isnumeric(x) && isscalar(x));
%    parser.addParameter('VideoWriterOptions', {});
    parser.parse(varargin{:});
    options = parser.Results;

    % set text interpreter
    set(groot, 'defaultTextInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');

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
    
    title('Steepest descent optimization process');
    xlabel('x1');
    ylabel('x2');
    if length(optimizationData.InitialData.OptimumCandidate) == 3
        zlabel('x3');
    end
    
    % Start video/gif
    if options.CreateVideoGifIterations
        filename = [options.SaveFolder, 'IterationProgression'];
        videoHandle = VideoWriter(filename, options.VideoWriterProfile);
        for i = 1:2:length(videoWriterOptions)
            videoHandle.(videoWriterOptions{i}) = videoWriterOptions{i+1};
        end
        open(videoHandle);
    end
    
    % Initialize graphic handles
    pointPrev = []; pointNext = [];
    arrowGrad = []; arrowSearch = [];
    stepLine = [];
    
    fixedArrowSize = options.ArrowScale;
    nIter = length(optimizationData.IterationData);
    
    for i = 1:nIter
        % Get design points and directions
        currentDesignVariables = optimizationData.IterationData(i).OptimumCandidatePrevious;
        nextDesignVariables = optimizationData.IterationData(i).OptimumCandidate;
        gradient = optimizationData.IterationData(i).Gradient;
        searchDirection = optimizationData.IterationData(i).SearchDirection;
    
        % Delete old elements if they exist and are valid
        if ~isempty(pointPrev) && isgraphics(pointPrev), delete(pointPrev); end
        if ~isempty(pointNext) && isgraphics(pointNext), delete(pointNext); end
        if ~isempty(arrowGrad) && isgraphics(arrowGrad), delete(arrowGrad); end
        if ~isempty(arrowSearch) && isgraphics(arrowSearch), delete(arrowSearch); end
        if ~isempty(stepLine) && isgraphics(stepLine), delete(stepLine); end
    
        % 2D plotting
        if length(currentDesignVariables) == 2
            pointPrev = plot(currentDesignVariables(1), currentDesignVariables(2), 'ko', ...
                             'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', 'Previous point');
            arrowGrad = quiver(currentDesignVariables(1), currentDesignVariables(2), ...
                               gradient(1)*fixedArrowSize, gradient(2)*fixedArrowSize, ...
                               0, 'Color', 'r', 'LineWidth', 1.5, 'MaxHeadSize', 2, 'DisplayName', 'Gradient');
            arrowSearch = quiver(currentDesignVariables(1), currentDesignVariables(2), ...
                                 searchDirection(1)*fixedArrowSize, searchDirection(2)*fixedArrowSize, ...
                                 0, 'Color', [0 0 0.5], 'LineWidth', 1.5, 'MaxHeadSize', 2, 'DisplayName', 'Search direction');
            pointNext = plot(nextDesignVariables(1), nextDesignVariables(2), 'go', ...
                             'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Current point');
            stepLine = plot([currentDesignVariables(1), nextDesignVariables(1)], ...
                            [currentDesignVariables(2), nextDesignVariables(2)], 'm--', ...
                            'LineWidth', 2, 'HandleVisibility', 'off');  % Hide from legend
        elseif length(currentDesignVariables) == 3
            pointPrev = plot3(currentDesignVariables(1), currentDesignVariables(2), currentDesignVariables(3), ...
                              'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
            arrowGrad = quiver3(currentDesignVariables(1), currentDesignVariables(2), currentDesignVariables(3), ...
                                gradient(1)*fixedArrowSize, gradient(2)*fixedArrowSize, gradient(3)*fixedArrowSize, ...
                                0, 'Color', 'r', 'LineWidth', 1.5, 'MaxHeadSize', 2);
            arrowSearch = quiver3(currentDesignVariables(1), currentDesignVariables(2), currentDesignVariables(3), ...
                                  searchDirection(1)*fixedArrowSize, searchDirection(2)*fixedArrowSize, searchDirection(3)*fixedArrowSize, ...
                                  0, 'Color', [0 0 0.5], 'LineWidth', 1.5, 'MaxHeadSize', 2);
            pointNext = plot3(nextDesignVariables(1), nextDesignVariables(2), nextDesignVariables(3), ...
                              'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
            stepLine = plot3([currentDesignVariables(1), nextDesignVariables(1)], ...
                             [currentDesignVariables(2), nextDesignVariables(2)], ...
                             [currentDesignVariables(3), nextDesignVariables(3)], 'm--', 'LineWidth', 2);
        end
    
        % Update axes
        set(gca, 'FontSize', 12, 'LineWidth', 1.5);
        legend('Location', 'northeast');
        xlim(xLimits); ylim(yLimits);
        if length(currentDesignVariables) == 3
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
    
    if options.CreateVideoGifIterations
        close(videoHandle);
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

    if isequal(options.LineSearchFunction, @line_search_backtracking)

        % Iterate over each optimization step
        for i = 1:nIter
            % Extract current design and gradient
            currentDesign = optimizationData.IterationData(i).OptimumCandidatePrevious;
            searchDirection = optimizationData.IterationData(i).SearchDirection;

            % Define step-size sweep
            stepsize_sweep = linspace(0, 1, 100);
            objective = zeros(1, length(stepsize_sweep));

            % Compute objective function values
            for j = 1:length(stepsize_sweep)
                objective(j) = evaluate_optimization_objective(objectiveFunction, currentDesign + stepsize_sweep(j) * searchDirection);
            end

            % Get step sizes used in the iteration (assuming stored in iteration data)
            stepsize_history = optimizationData.IterationData(i).StepSizeData; % Extracted step sizes
            final_stepsize = stepsize_history(end); % Last step size is α*

            % Clear figure and plot new iteration
            clf; hold on;

            % Plot the step-size sweep vs objective function
            plot(stepsize_sweep, objective, 'b-', 'LineWidth', 2);
            xlabel('Step Size (\alpha)');
            ylabel('Objective Function Value');
            title(['Iteration ', num2str(i)]);
            grid on;

            % Initialize filename for this iteration's GIF
            gif_filename = fullfile(lineSearchFolder, sprintf('LineSearch_Iteration_%d.gif', i));

            % Plot step-size intersections with different colors
            for j = 1:length(stepsize_history)
                objectiveVariable = evaluate_optimization_objective(objectiveFunction, currentDesign + stepsize_history(j) * searchDirection);

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
                plot(stepsize_history(j), objectiveVariable, 'o', ...
                     'MarkerSize', markerSize, 'MarkerEdgeColor', 'k', ...
                     'MarkerFaceColor', markerColor, 'LineWidth', 1.5);

                % Add xline as a reference
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
            currentDesign = optimizationData.IterationData(i).OptimumCandidatePrevious;
            searchDirection = optimizationData.IterationData(i).SearchDirection;
    
            % Define step-size sweep
            stepsize_sweep = linspace(0, 3, 100); % Sweeping from 0 to 1
            objective = zeros(1, length(stepsize_sweep));
    
            % Compute objective function values
            for j = 1:length(stepsize_sweep)
                objective(j) = evaluate_optimization_objective(objectiveFunction, currentDesign + stepsize_sweep(j) * searchDirection);
            end
    
            % Get step sizes used in the iteration
            stepsize_history = optimizationData.IterationData(i).StepSizeData; % Extracted step sizes
            final_stepsize = stepsize_history(end); % Last step size is α*
    
            % Plot the line search (Objective function curve)
            clf; hold on;
            plot(stepsize_sweep, objective, 'b-', 'LineWidth', 2); % Plot objective curve
            xlabel('Step size (\alpha)');
            ylabel('Objective function value');
            title(['Iteration ', num2str(i)]);
            grid on;
    
            % Initialize filename for this iteration's GIF
            gif_filename = fullfile(lineSearchFolder, sprintf('LineSearch_Iteration_%d.gif', i));
    
            % Plot step-size intersections with different colors
            for j = 1:length(stepsize_history)
                % Find corresponding objective function value
                objectiveVariable = evaluate_optimization_objective(objectiveFunction, currentDesign + stepsize_history(j) * searchDirection);
    
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
                plot(stepsize_history(j), objectiveVariable, 'o', ...
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
    set(gcf, 'Color', 'w');  % Set the figure background to white
    hold on;
    grid on;
    axis equal;
    
    % Extract problem bounds
    xLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(1), optimizationData.ProblemData.DesignSpaceUpperBound(1)];
    yLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(2), optimizationData.ProblemData.DesignSpaceUpperBound(2)];
    
    % Plot contour once (for 2D only)
    if length(optimizationData.InitialData.OptimumCandidate) == 2
        hContour = plot_point_optimization_contour_2d(figureHandleOptimization, optimizationData.ProblemData(1).ObjectiveFunction, contourDesignSpace, ...
            'DesignVariableNames', options.DesignVariableNames, ...
            'AxisFontSize', options.AxisFontSize, ...
            options.PlotContour2DOptions{:});

        % Hide all graphics objects created in this function from legend
        if isgraphics(hContour)
            set(findall(hContour, '-property', 'HandleVisibility'), 'HandleVisibility', 'off');
        end
    end
    
    % Plot the contour lines
    xlabel('x1');
    ylabel('x2');
    title('Quasi newton optimization progress');
    
    % Get number of iterations
    nIter = length(optimizationData.IterationData);
    
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
    imwrite(A, map, gif_filename, 'gif', 'LoopCount', inf, 'DelayTime', frameDelay, 'WriteMode', 'overwrite');
    
    % Pause before next iteration
    pause(frameDelay);
    
    % Loop through iterations, displaying points one by one
    for i = 1:nIter-1
        % Extract current and next points
        currentPoint = optimizationData.IterationData(i).OptimumCandidatePrevious;
        nextPoint = optimizationData.IterationData(i).OptimumCandidate;
        
        % First, plot the current point (red)
        plot(currentPoint(1), currentPoint(2), 'ko', ...
             'MarkerSize', 10, 'MarkerFaceColor', 'k', 'DisplayName', 'Current point'); % black for current point
        
        % Store previous points for connection
        previousX = [previousX, currentPoint(1)];
        previousY = [previousY, currentPoint(2)];
        
        % Pause to make the appearance gradual
        pause(frameDelay);
        
        % After the pause, plot the connection line (magenta) between the current and next point
        plot([previousX(end), nextPoint(1)], [previousY(end), nextPoint(2)], 'm--', 'LineWidth', 2); % Magenta line
        
        % Capture the frame for GIF
        frame = getframe(gcf);
        img = frame2im(frame);
        [A, map] = rgb2ind(img, 256);
    
        % Save as GIF (subsequent frames)
        imwrite(A, map, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frameDelay);
        
        % Pause before next iteration
        pause(frameDelay);
    end
    
    % Plot the final optimum point in blue
    finalPoint = optimizationData.IterationData(nIter).OptimumCandidate;
    plot(finalPoint(1), finalPoint(2), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b'); % Blue for final optimum
    
    % Connect last iteration to final point
    plot([previousX(end), finalPoint(1)], [previousY(end), finalPoint(2)], 'm--', 'LineWidth', 2); % magenta connecting line
    
    % Capture the final frame for GIF
    frame = getframe(gcf);
    img = frame2im(frame);
    [A, map] = rgb2ind(img, 256);
    
    % Save final frame in the GIF
    imwrite(A, map, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frameDelay);
    
    % Close video handle
    if (options.CreateVideoGifIterations)
        close(videoHandle);
    end
end




