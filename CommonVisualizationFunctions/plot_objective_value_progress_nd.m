function plot_objective_value_progress_nd(optimizationData, varargin)
%PLOT_OBJECTIVE_VALUE_PROGRESS_ND Plot objective value progression in n-dimensional design space.
%
%   plot_objective_value_progress_nd(optimizationData, ...)
%   plots the objective value progression in an n-dimensional design space.
%
%   INPUTS:
%       optimizationData : struct containing optimization data
%
%   NAME-VALUE PAIR ARGUMENTS:
%       Algorithm : string
%           Algorithm used for optimization.
%       SaveFolder : string
%           Folder to save the plot.
%       SaveFigureOptions : cell array
%           Options for saving the figure.
%       CreateVideoGifIterations : logical
%           Whether to create a video of the iterations.
%       FramesPerSecond : double
%           Frames per second for the video.
%       ImwriteOptions : cell array
%           Options for imwrite.
%       VideoWriterProfile : string
%           Profile for the video writer.
%       VideoWriterOptions : cell array
%           Options for the video writer.
%       ShowContour : logical
%           Whether to show the contour plot.
%       ContourLevels : double
%           Number of contour levels to show.
%       ContourOptions : cell array
%           Options for the contour plot.
%       IncludeLegend : logical
%           Whether to include the legend.
%       IncludeArrow : logical
%           Whether to include the arrow.
%       LegendOptions : cell array
%           Options for the legend.
%       LineWidth : double
%           Width of the line.
%       Color : string
%           Color of the line.
%       MarkerEdgeColor : string
%           Color of the marker edge.
%       Marker : string
%           Marker of the line.
%       AxisFontSize : double
%           Font size of the axis.
%       TitleFontSize : double
%           Font size of the title.
%       LineStyle : string
%           Style of the line.
%       LineType : string
%           Type of the line.
%       ShowGrid : string
%           Whether to show the grid.
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
    p = inputParser;
    addParameter(p, 'Algorithm', 'Particle Swarm');
    addParameter(p, 'SaveFolder', []);
    addParameter(p, 'SaveFigureOptions', {});
    addParameter(p, 'CreateVideoGifIterations', true);
    addParameter(p, 'FramesPerSecond', 5);
    addParameter(p, 'ImwriteOptions', {});
    addParameter(p, 'VideoWriterProfile', 'Motion JPEG 2000');
    addParameter(p, 'VideoWriterOptions', {'LosslessCompression', true});
    addParameter(p, 'ShowContour', true);
    addParameter(p, 'ContourLevels', 10);
    addParameter(p, 'ContourOptions', {});
    addParameter(p, 'IncludeLegend', true);
    addParameter(p, 'IncludeArrow', true);
    addParameter(p, 'LegendOptions', {'location', 'southeast', 'FontSize', 6});
    addParameter(p, 'LineWidth', 0.5);
    addParameter(p, 'Color', 'b');
    addParameter(p, 'MarkerEdgeColor', 'b');
    addParameter(p, 'Marker', 'o');
    addParameter(p, 'AxisFontSize', 12);
    addParameter(p, 'TitleFontSize', 14);
    addParameter(p, 'LineStyle', '-');
    addParameter(p, 'LineType', 'line');
    addParameter(p, 'ShowGrid', 'on');
    parse(p, varargin{:});
    options = p.Results;

    defaultVideoWriterOptions = {'FrameRate',options.FramesPerSecond};
    [~,videoWriterOptions] = merge_name_value_pair_argument(defaultVideoWriterOptions,options.VideoWriterOptions);
    defaultImwriteOptions = {};
    [~,imwriteOptions] = merge_name_value_pair_argument(defaultImwriteOptions,options.ImwriteOptions);
    frameDelay = 1/options.FramesPerSecond;
    defaultLegendOptions = {'location','southeast','FontSize', 6};
    [~,legendOptions] = merge_name_value_pair_argument(defaultLegendOptions,options.LegendOptions);
    
    % Get dimensions and pairs
    numberDimensions = length(optimizationData.ProblemData.DesignSpaceLowerBound);
    dimensionPairs = nchoosek(1:numberDimensions, 2);
    numberPairs = size(dimensionPairs, 1);
    
    % Optimization iterations
    nIter = length(optimizationData.IterationData);
    if nIter == 0
        warning('No iteration data available. Exiting.');
        return;
    end
    
    % Gather optimum candidates over iterations (including initial)
    optimumCandidate = [optimizationData.InitialData.OptimumCandidate; vertcat(optimizationData.IterationData(:).OptimumCandidate)];
    optimumCandidateObjectiveValue = [optimizationData.InitialData.OptimumCandidateObjectiveValue; vertcat(optimizationData.IterationData(:).OptimumCandidateObjectiveValue)];
    
    gifFilenames = cell(numberPairs,1);
    % Prepare video writers if saving
    if options.CreateVideoGifIterations && ~isempty(options.SaveFolder)
        videoHandles = cell(numberPairs, 1);
        for iPair = 1:numberPairs
            filename = fullfile(options.SaveFolder, sprintf('ObjectiveValueProgression_x%d_x%d.avi', dimensionPairs(iPair, 1), dimensionPairs(iPair, 2)));
            videoHandle = VideoWriter(filename, options.VideoWriterProfile);
            for idx = 1:2:length(videoWriterOptions)
                try
                    videoHandle.(videoWriterOptions{i}) = videoWriterOptions{i+1};
                catch
                    warning('Could not set VideoWriter property');
                end
            end
            open(videoHandle);
            videoHandles{iPair} = videoHandle;
            gifFilenames{iPair} = [filename '.gif'];
        end
    else
        videoHandles = [];
    end
    
    % Precompute contours if needed
    if options.ShowContour
        contourData = cell(numberPairs, 1);
        for iPair = 1:numberPairs
            dimension1 = dimensionPairs(iPair, 1);
            dimension2 = dimensionPairs(iPair, 2);
            
            xLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(dimension1), optimizationData.ProblemData.DesignSpaceUpperBound(dimension1)];
            yLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(dimension2), optimizationData.ProblemData.DesignSpaceUpperBound(dimension2)];
            
            [x,y] = meshgrid(linspace(xLimits(1), xLimits(2), 100), linspace(yLimits(1), yLimits(2), 100));
            
            % Build full input points with dims dim1 and dim2 replaced by xg, yg grid
            nPoints = numel(x);
            inputGrid = repmat(optimumCandidate(1, :), nPoints, 1);
            inputGrid(:, dimension1) = x(:);
            inputGrid(:, dimension2) = y(:);
            
            zVals = evaluate_optimization_objective(optimizationData.ProblemData.ObjectiveFunction, inputGrid);
            zVals = reshape(zVals, size(x));
            
            contourData{iPair}.x = x;
            contourData{iPair}.y = y;
            contourData{iPair}.z = zVals;
            contourData{iPair}.cmin = min(zVals(:));
            contourData{iPair}.cmax = max(zVals(:));
        end
    end
    
    % -------------------------------------------------
    % Creating layout for n dimenison
    % -------------------------------------------------
    
    % Number of columns and rows for scatter plots
    nColumns = 2;
    nRows = ceil(numberPairs / nColumns);
    
    % Desired gaps (normalized units)
    gapBetweenScatterPlotsX = 0.1; % horizontal gap between scatter plots
    gapBetweenScatterPlotsY = 0.1; % vertical gap between scatter plots
    gapBetweenScatterAndLine = 0.4; % gap between 2nd scatter column and line plot
    
    % Total width for scatter plots area (before line plot)
    scatterAreaWidth = 1 - 0.05 - gapBetweenScatterAndLine - 0.05; % 0.05 left margin + 0.05 right margin approx
    scatterAreaHeight = 0.85; % approx height for scatter plots area
    
    % Compute width and height for each scatter plot
    scatterWidth = (scatterAreaWidth - (nColumns-1)*gapBetweenScatterPlotsX) / nColumns;
    scatterHeight = (scatterAreaHeight - (nRows-1)*gapBetweenScatterPlotsY) / nRows;
    
    % Left margin from figure left
    leftMargin = 0.05;
    
    % --- Create figure
    figureHandle = figure('Color', 'w', 'MenuBar', 'none', 'ToolBar', 'none');
    set(figureHandle, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    
    % -------------------------------------------------
    % Initializing Line plot
    % -------------------------------------------------
    lineAxesPosition = [0.7 0.1 0.2 0.8]; % width 0.3, height 0.8 as example
    lineAxes = axes('Position', lineAxesPosition);
    hold(lineAxes, 'on');
    xlabel(lineAxes, 'Iteration (k)', 'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
    ylabel(lineAxes, 'Optimal Objective Value ($f(\mathbf{x})$)', 'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
    grid(lineAxes, options.ShowGrid);
    xlim(lineAxes, [1 nIter]);
    ylim(lineAxes, [min(optimumCandidateObjectiveValue) max(optimumCandidateObjectiveValue)]);
    title(lineAxes, 'Objective Value Evolution', 'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
    handleObjectiveValue = plot(lineAxes, NaN, NaN, ...
        'LineWidth', options.LineWidth, 'Color', options.Color, ...
        'Marker', options.Marker, 'LineStyle', options.LineStyle);
    
    % -------------------------------------------------
    % Initializing Scatter plot
    % -------------------------------------------------
    scatterAxes = gobjects(numberPairs, 1);
    for iPair = 1:numberPairs
        if numberPairs == 1
                scatterAxes(iPair) = axes('Position', [0.1 0.15 0.5 0.5]); % Enlarged position
        else
            % Calculate row and column indices
            row = ceil(iPair / nColumns);
            col = mod(iPair-1, nColumns) + 1;
            
            % Compute position for each scatter plot
            leftPos = leftMargin + (col-1)*(scatterWidth + gapBetweenScatterPlotsX);
            % Bottom positions count from bottom up (0 is bottom), so invert row index
            bottomPos = 1 - (row * scatterHeight + (row - 1)*gapBetweenScatterPlotsY) - 0.1; % 0.1 bottom margin approx
            
            % Create axes
            scatterAxes(iPair) = axes('Position', [leftPos bottomPos scatterWidth scatterHeight]);
        end
        hold(scatterAxes(iPair), 'on');

        dimension1 = dimensionPairs(iPair, 1);
        dimension2 = dimensionPairs(iPair, 2);
        
        % Plot contour
        if options.ShowContour
            [~, h] = contourf(scatterAxes(iPair), contourData{iPair}.x, contourData{iPair}.y, contourData{iPair}.z, ...
                options.ContourLevels, options.ContourOptions{:});
            contourHandles(iPair) = h;
            colormap(scatterAxes(iPair), flipud(gray(256)));
            caxis(scatterAxes(iPair), [contourData{iPair}.cmin contourData{iPair}.cmax]);
            colorbar(scatterAxes(iPair));
        end
        
        % Initialize scatter plots with empty data
        scatterPopHandles(iPair) = scatter(scatterAxes(iPair), NaN, NaN, 36, 'b', 'filled', 'MarkerEdgeColor', options.MarkerEdgeColor);
        scatterOptHandles(iPair) = scatter(scatterAxes(iPair), NaN, NaN, 100, 'r', 'filled', 'Marker', 'p', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        
        xlabel(scatterAxes(iPair), sprintf('x%d', dimension1), 'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
        ylabel(scatterAxes(iPair), sprintf('x%d', dimension2), 'Interpreter', 'latex', 'FontSize', options.AxisFontSize);
        grid(scatterAxes(iPair), options.ShowGrid);
        axis(scatterAxes(iPair), [optimizationData.ProblemData.DesignSpaceLowerBound(dimension1) optimizationData.ProblemData.DesignSpaceUpperBound(dimension1) ...
            optimizationData.ProblemData.DesignSpaceLowerBound(dimension2) optimizationData.ProblemData.DesignSpaceUpperBound(dimension2)]);
        title(scatterAxes(iPair), sprintf('Current Population: x%d vs x%d', dimension1, dimension2), 'Interpreter', 'latex', 'FontSize', options.TitleFontSize);
        axis equal;
        if options.IncludeLegend
            lgd = legend(scatterAxes(iPair), [scatterPopHandles(iPair), scatterOptHandles(iPair)], ...
                         {'Population', 'Optimal Candidate'}, legendOptions{:});
            set(lgd, 'Box', 'on', 'Color', [0.9 0.9 0.9], 'EdgeColor', 'k');
        end
    end
    
    % -------------------------------------------------
    % Update loop per iteration
    % -------------------------------------------------

    for optimizationIteration = 1:nIter
        population = optimizationData.IterationData(optimizationIteration).PopulationDesignNew;
        
        % Update scatter points for all pairs
        for iPair = 1:numberPairs
            dimension1 = dimensionPairs(iPair, 1);
            dimension2 = dimensionPairs(iPair, 2);
            set(scatterPopHandles(iPair), 'XData', population(:, dimension1), 'YData', population(:, dimension2));
            set(scatterOptHandles(iPair), 'XData', optimumCandidate(optimizationIteration+1, dimension1), 'YData', optimumCandidate(optimizationIteration+1, dimension2));
        end

        sgtitle(sprintf('Algorithm: %s Iteration %d, Population Size %d', ...
            options.Algorithm, optimizationIteration, optimizationData.ProblemData.Options.PopulationSize)); % Adds overarching title
        
        % Update line plot
        set(handleObjectiveValue, 'XData', 1:optimizationIteration, 'YData', optimumCandidateObjectiveValue(2:optimizationIteration+1));
        
        drawnow;
        
        if options.CreateVideoGifIterations && ~isempty(options.SaveFolder)
            currentFrame = getframe(gcf);
            [A, map] = rgb2ind(frame2im(currentFrame), 256);
        
            for iPair = 1:numberPairs
                % Resize frame if not 560x420
                frameImage = frame2im(currentFrame);
                if size(frameImage, 1) ~= 420 || size(frameImage, 2) ~= 560
                    frameImage = imresize(frameImage, [420, 560]);
                    currentFrame = im2frame(frameImage);
                end
        
                writeVideo(videoHandles{iPair}, currentFrame);
        
                if optimizationIteration == 1
                    imwrite(A, map, gifFilenames{iPair}, 'gif', ...
                        'LoopCount', inf, ...
                        'DelayTime', frameDelay, ...
                        imwriteOptions{:});
                else
                    imwrite(A, map, gifFilenames{iPair}, 'gif', ...
                        'WriteMode', 'append', ...
                        'DelayTime', frameDelay, ...
                        imwriteOptions{:});
                end
            end
        end
    end
    
    % Close videos if any
    if options.CreateVideoGifIterations && ~isempty(options.SaveFolder)
        for iPair = 1:numberPairs
            close(videoHandles{iPair});
        end
    end
end
