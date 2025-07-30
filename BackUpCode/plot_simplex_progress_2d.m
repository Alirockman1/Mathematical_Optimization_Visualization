function plot_simplex_progress_2d(optimizationData, objectiveFunction, varargin)
    parser = inputParser;
    parser.addParameter('SaveFolder', []);
    parser.addParameter('SaveFigureOptions', {});
    parser.addParameter('CreateVideoGifIterations', true);
    parser.addParameter('FramesPerSecond', 5);
    parser.addParameter('ImwriteOptions', {});
    parser.addParameter('VideoWriterProfile', 'Motion JPEG 2000');
    parser.addParameter('VideoWriterOptions', {'LosslessCompression', true});
    parser.addParameter('AxisFontSize', 12, @(x) isnumeric(x) && isscalar(x));
    parser.addParameter('DesignVariableNames',{});
    parser.addParameter('UseDesignSpaceAsAxisLimits',true);
    parser.addParameter('PlotContour2DUseDesignSpaceLimit',true);
    parser.addParameter('PlotContour2DOptions',{});
    parser.addParameter('Title','Linear programming (simplex) optimization progress',@ischar);
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

    % Plot contour once (for 2D only)
    if length(optimizationData.InitialData.OptimumCandidate) == 2
        if options.PlotContour2DUseDesignSpaceLimit
            contourDesignSpace = [optimizationData.ProblemData(1).DesignSpaceLowerBound(:)';
                optimizationData.ProblemData(1).DesignSpaceUpperBound(:)'];
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

    % Plot the contour lines
    xlabel('x1');
    ylabel('x2');
    title(options.Title);
    if length(optimizationData.InitialData.OptimumCandidate) == 3
        zlabel('x3');
    end

    % Get number of iterations
    nIter = length(optimizationData.IterationData);

    % Initialization
    previousX = [];
    previousY = [];

    % Start video/gif
    if options.CreateVideoGifIterations
        filename = [options.SaveFolder, 'IterationProgression'];
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




