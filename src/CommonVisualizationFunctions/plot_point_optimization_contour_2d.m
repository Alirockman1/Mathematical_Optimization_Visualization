function plotHandle = plot_point_optimization_contour_2d(figureHandle, objectiveFunction, designSpace, varargin)
%PLOT_POINT_OPTIMIZATION_CONTOUR_2D Generate 2D contour plot of optimization objective function.
%
%   PLOTHANDLE = PLOT_POINT_OPTIMIZATION_CONTOUR_2D(FIGUREHANDLE, OBJECTIVEFUNCTION, ...
%                   DESIGNSPACE, ...) creates a contour plot of a 2D optimization
%                   landscape within specified design space bounds.
%
%   INPUTS:
%       FIGUREHANDLE       - Handle to target figure for plotting
%       OBJECTIVEFUNCTION  - Function handle to objective function (must accept N×2 matrix)
%       DESIGNSPACE        - 2×2 matrix specifying bounds:
%           [lb_x, lb_y;   - Lower bounds for design variables
%            ub_x, ub_y]   - Upper bounds for design variables
%
%   NAME-VALUE PAIR ARGUMENTS:
%       'GridRelativeStepSize' - Grid resolution as fraction of design space (default: 0.01)
%       'ContourLevels'       - Number or specific values for contour lines (default: 10)
%       'ContourOptions'      - Cell array of additional contour plot options (default: {})
%       'DesignVariableNames' - Custom names for design variables (default: x₁,x₂)
%       'AxisFontSize'        - Font size for axis labels (default: 14)
%
%   OUTPUTS:
%       PLOTHANDLE          - Handle to contour plot object (optional)
%
%   FUNCTIONALITY:
%       1. Creates uniform grid across specified design space
%       2. Evaluates objective function at grid points
%       3. Generates contour plot with automatic or specified levels
%       4. Applies LaTeX formatting to axis labels
%       5. Adds colormap and colorbar for objective values
%
%   NOTES:
%       - Uses jet colormap by default
%       - Automatically calculates 10 contour levels if not specified
%       - Maintains aspect ratio of design space
%       - Preserves existing figure content when using hold on
%
%   DEPENDENCIES:
%       - evaluate_optimization_objective.m
%       - merge_name_value_pair_argument.m
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor)
%   Copyright 2025 Ali Abbas Kapadia (Author) 
%   Copyright 2025 Gaurav Vaibhav (Author)
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
    parser.addParameter('GridRelativeStepSize',0.01);
    parser.addParameter('ContourLevels',[]);
    parser.addParameter('ContourOptions',{});
    parser.addParameter('DesignVariableNames',{});
    parser.addParameter('AxisFontSize', 14, @(x) isnumeric(x) && isscalar(x));
    parser.parse(varargin{:});
    options = parser.Results;

    % design variable names
    if(isempty(options.DesignVariableNames))
        for i=1:length(designSpace(1,:))
            options.DesignVariableNames{i} = sprintf('$$x_{%d}$$',i);
        end
    end
 
    % make the intervals for each variable: interval = lb + grid*(ub - lb)
    xInterval = designSpace(1,1) + (0:options.GridRelativeStepSize:1)*(designSpace(2,1)-designSpace(1,1));
    yInterval = designSpace(1,2) + (0:options.GridRelativeStepSize:1)*(designSpace(2,2)-designSpace(1,2));
 
    % generate grid for predictions at finer sample rate
    [xGrid, yGrid] = meshgrid(xInterval,yInterval);
    fullGrid = [xGrid(:),yGrid(:)];
    [objectiveValue,~] = evaluate_optimization_objective(objectiveFunction,fullGrid);
    %objectiveValue = objectiveFunction(fullGrid);

    % reshape to same grid size as the input
    objectiveValue = reshape(objectiveValue, size(xGrid));

    % get desired level of contours
    if(isempty(options.ContourLevels))
        levels = 10;
    elseif(isscalar(options.ContourLevels))
        levels = linspace(min(objectiveValue(:)), max(objectiveValue(:)), options.ContourLevels);
    else
        levels = options.ContourLevels;
    end

    %plot options
    defaultPlotOptions = {};
    [~,contourOptions] = merge_name_value_pair_argument(defaultPlotOptions,options.ContourOptions);

    % plot decision boundary surface
    figure(figureHandle);
    hold all;
    [~,plotHandle] = contour(xGrid,yGrid,objectiveValue,levels,contourOptions{:});
    xlabel(options.DesignVariableNames{1},'interpreter','latex','FontSize',options.AxisFontSize);
    ylabel(options.DesignVariableNames{2},'interpreter','latex','FontSize',options.AxisFontSize);
    colormap('jet');
    colorbar; % Optionally add a color bar
    grid on;

    if(nargout<1)
        clear plotHandle;
    end
end
