function plot_objective_value_progress_2d(optimizationData,varargin)
% PLOT_OBJECTIVE_VALUE_PROGRESS_2D Visualizes 2D optimization progress with objective value evolution
%
%   PLOT_OBJECTIVE_VALUE_PROGRESS_2D(OPTIMIZATIONDATA) creates a dual-panel 
%   visualization showing:
%     - Left panel: Population distribution in design space with contour plot
%     - Right panel: Objective value convergence over iterations
%     - Optional animation output as GIF/video
%
%   INPUTS:
%       optimizationData - Structure containing optimization results with:
%           .ProblemData - Problem definition including bounds and objective
%           .IterationData - Array of optimization history
%           .InitialData - Initial population and values
%
%   OPTIONAL PARAMETERS (name-value pairs):
%       'Algorithm' - Name of optimization algorithm (default: 'Differential Algorithm')
%       'SaveFolder' - Path to save output figures/videos (default: [])
%       'SaveFigureOptions' - Options for figure saving (cell array)
%       'CreateVideoGifIterations' - Enable GIF/video creation (default: true)
%       'FramesPerSecond' - Animation frame rate (default: 5)
%       'ImwriteOptions' - Options for GIF creation (cell array)
%       'VideoWriterProfile' - Video profile (default: 'Motion JPEG 2000')
%       'VideoWriterOptions' - Video writing options (cell array)
%       'ShowContour' - Show objective function contours (default: true)
%       'ContourLevels' - Number of contour levels (default: 10)
%       'ContourOptions' - Additional contour options (cell array)
%       'IncludeLegend' - Show legend (default: true)
%       'IncludeArrow' - Show direction arrows (default: true)
%       Visual style parameters:
%           'LineWidth', 'Color', 'MarkerEdgeColor', 'Marker', 
%           'AxisFontSize', 'LineStyle', 'LineType', 'ShowGrid'
%       'LegendOptions' - Legend customization options (cell array)
%
%   OUTPUTS:
%       Figure with:
%       1. Left panel: 
%          - Population distribution (blue points)
%          - Current optimum candidate (red point)
%          - Optional objective function contour
%       2. Right panel:
%          - Objective value convergence curve
%       3. Optionally saves animation as GIF/video
%
%   KEY FEATURES:
%       - Dual-panel view showing both design space and convergence
%       - Interactive animation through iterations
%       - Professional formatting with LaTeX labels
%       - Flexible visualization customization
%       - Publication-quality output options
%
%   ALGORITHM VISUALIZATION:
%       The visualization helps understand:
%       1. How the population explores the design space
%       2. Relationship between design space location and objective value
%       3. Convergence behavior of the optimization
%
%   SEE ALSO:
%       plot_differential_evolution_progress_nd_slider, 
%       plot_contour_objective_constraint_results_slider
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
    parser.addParameter('Algorithm','Differential Algorithm', @(x) ischar(x) || isstring(x));
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
    parser.addParameter('LineWidth', 0.5, @(x) isnumeric(x) && isscalar(x));
    parser.addParameter('Color', 'b', @(x) ischar(x) || isstring(x));
    parser.addParameter('MarkerEdgeColor', 'b', @(x) ischar(x) || isstring(x));
    parser.addParameter('Marker', 'o', @(x) ischar(x) || isstring(x));
    parser.addParameter('AxisFontSize', 12, @(x) isnumeric(x) && isscalar(x));
    parser.addParameter('LineStyle', '-', @(x) ischar(x) || isstring(x));
    parser.addParameter('LineType', 'line', @(x) ischar(x) || isstring(x));
    parser.addParameter('ShowGrid', 'on', @(x) islogical(x) || ischar(x) && any(strcmpi(x, {'on', 'off'})));
    parser.parse(varargin{:});
    options = parser.Results;

    defaultVideoWriterOptions = {'FrameRate',options.FramesPerSecond};
    [~,videoWriterOptions] = merge_name_value_pair_argument(defaultVideoWriterOptions,options.VideoWriterOptions);
    defaultImwriteOptions = {};
    [~,imwriteOptions] = merge_name_value_pair_argument(defaultImwriteOptions,options.ImwriteOptions);
    frameDelay = 1/options.FramesPerSecond;
    defaultLegendOptions = {'location','southeast','FontSize', 6};
    [~,legendOptions] = merge_name_value_pair_argument(defaultLegendOptions,options.LegendOptions);
    
    % number of iterations
	nIter = length(optimizationData.IterationData);

    optimumCandidate = [optimizationData.InitialData.OptimumCandidate;...
        vertcat(optimizationData.IterationData(:).OptimumCandidate)];
    optimumCandidateObjectiveValue = [optimizationData.InitialData.OptimumCandidateObjectiveValue;...
        vertcat(optimizationData.IterationData(:).OptimumCandidateObjectiveValue)];


	% start video/gif
	if(options.CreateVideoGifIterations)
		filename = [options.SaveFolder,'ObjectiveValueProgression'];
		videoHandle = VideoWriter(filename,options.VideoWriterProfile);
	    for i=1:2:length(videoWriterOptions)
	        videoHandle.(videoWriterOptions{i}) = videoWriterOptions{i+1};
	    end
	    open(videoHandle);
	end

	figure;
	hold on;
	axis([optimizationData.ProblemData.DesignSpaceLowerBound(1) optimizationData.ProblemData.DesignSpaceUpperBound(1)...
		optimizationData.ProblemData.DesignSpaceLowerBound(2) optimizationData.ProblemData.DesignSpaceUpperBound(2)]);
    xLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(1), optimizationData.ProblemData.DesignSpaceUpperBound(1)];
    yLimits = [optimizationData.ProblemData.DesignSpaceLowerBound(2), optimizationData.ProblemData.DesignSpaceUpperBound(2)];

    if options.ShowContour
        [x,y] = meshgrid(linspace(xLimits(1),xLimits(2),100),linspace(yLimits(1),yLimits(2),100));
        fullGrid = [x(:),y(:)];
        z = evaluate_optimization_objective(optimizationData.ProblemData.ObjectiveFunction,fullGrid);
        z = reshape(z, size(x));  
    end

	grid('minor');
	set(gcf, 'MenuBar', 'none');
    set(gcf, 'ToolBar', 'none');
    set(gcf, 'Color', 'w');  

	for i=1:nIter
		if(i==1)
			populationPrevious = optimizationData.InitialData.PopulationDesignInitial;
		else
			populationPrevious = optimizationData.IterationData(i-1).PopulationDesignNew;
        end
        % Add a single title for the entire figure
        handleIterationTitle = sgtitle(sprintf('%s Population Size %d, Iteration Index (k) = %d',options.Algorithm, optimizationData.ProblemData.Options.PopulationSize, i)); % Adds overarching title
        
        % Subframe 1: Current Population
        ax1 = subplot(1, 2, 1); % Top-left corner
        hold on;
        if options.ShowContour
            contourf(x,y,z,options.ContourLevels,options.ContourOptions{:});
            colormap(ax1,flipud(gray(256)));
            cmin = min(z(:));
            cmax = max(z(:));
            clim([cmin cmax]); 
            cbar = colorbar;
            cbar.Label.String = 'Objective Value';
        end
        size(optimumCandidate(i,:));
        size(populationPrevious);
        indexOptimum = ismember(populationPrevious,optimumCandidate(i,:),'rows');
        handleCurrentPopulation = plot(populationPrevious(~indexOptimum, 1), populationPrevious(~indexOptimum, 2),'.', 'Color','b','MarkerSize', 10);
        handleOptimumCandidate = plot(optimumCandidate(i,1),optimumCandidate(i,2),'.','Color','r','MarkerSize',12);
        title('Current Population','FontSize',options.AxisFontSize);
        xlabel('$$x_1$$','interpreter','latex','FontSize',options.AxisFontSize); 
        ylabel('$$x_2$$','interpreter','latex','FontSize',options.AxisFontSize);
        xlim(xLimits);
        ylim(yLimits);
        axis equal;
        grid(options.ShowGrid);
        if options.IncludeLegend
            legend(handleOptimumCandidate, {'Optimal Candidate'}, legendOptions{:});
            set(legend, 'Box', 'on', 'Color', [0.9 0.9 0.9], 'EdgeColor', 'k');
        end

        % Subframe 2: Mutation Step
        ax2 = subplot(1, 2, 2); % Top-right corner
        hold on;
        if i==1
            handleObjectiveValue = plot(i, optimumCandidateObjectiveValue(i), 'LineWidth', options.LineWidth, 'Color', options.Color, ...
                    'Marker', options.Marker, 'LineStyle', options.LineStyle);
        else
            handleObjectiveValue = plot([i-1;i], [optimumCandidateObjectiveValue(i-1);optimumCandidateObjectiveValue(i)],...
            'LineWidth', options.LineWidth, 'Color', options.Color,'Marker', options.Marker, 'LineStyle', options.LineStyle);
        end
        
        title('Optimal Objective Value Evolution','interpreter','latex','FontSize',options.AxisFontSize);
        xlabel('Iteration Index (k)','interpreter','latex','FontSize',options.AxisFontSize); 
        ylabel('Optimal Objective Value ($f(\mathbf{x})$)','interpreter','latex','FontSize',options.AxisFontSize);
        xlim([1 nIter]);
        ylim([min(optimumCandidateObjectiveValue) max(optimumCandidateObjectiveValue)]);
        grid(options.ShowGrid);
        
        % Adjust all subplots to have same color range
        set(ax1, 'CLim', [cmin cmax]);

    % Save video/GIF frame after all subframes for this iteration
		if(options.CreateVideoGifIterations)
			% create new video frame 
			currentFrame = getframe(gcf);
        	writeVideo(videoHandle,currentFrame);

			% create new gif frame
			[A,map] = rgb2ind(frame2im(getframe(gcf)),256);
	        if(i==1)
	            imwrite(A,map,...
	                [filename,'.gif'],'gif',...
	                'LoopCount',inf,...
	                'DelayTime',frameDelay,...
	                imwriteOptions{:});
	        else
	            imwrite(A,map,...
	                [filename,'.gif'],'gif',...
	                'WriteMode','append', ...
	                'DelayTime',frameDelay,...
	                imwriteOptions{:});
	        end
        end

		delete(handleCurrentPopulation);
        delete(handleOptimumCandidate);
        delete(handleIterationTitle);
	end
	% reset menu/tool bars
    set(gcf, 'MenuBar', 'figure');
    set(gcf, 'ToolBar', 'auto');

    if(options.CreateVideoGifIterations)
    	close(videoHandle);
    end
end