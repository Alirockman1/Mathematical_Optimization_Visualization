%setup_optimization_visualization_toolbox
%
%   Copyright 2025 Eduardo Rodrigues Della Noce (Supervisor, Author)
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


%% Cleanup
close all;
fclose all;
clear all;
clc;
more off;
diary off;


%% get folder of setup
setupFilePath = mfilename('fullpath');
if contains(setupFilePath,'LiveEditorEvaluationHelper')
    setupFilePath = matlab.desktop.editor.getActiveFilename;
end
toolboxFolderPath = replace(erase(setupFilePath,mfilename),'\','/');
addpath(toolboxFolderPath);
clear setupFilePath


%% Add Folders (and potentially Subfolders) to Path
% 'foldername',includeSubfolder
addDirectory = {...
    ... % public folders
    'CommonConvergenceFunctions/',true;...
    'CommonEvaluationFunctions/',true;...
    'CommonVisualizationFunctions/',true;...
    'ConjugateGradient/',true;...
    'DifferentialEvolution/',true;...
    'GeneticAlgorithm/',true;...
    'LinearPrograming/',true;...
    'SequentialLinearPrograming/',true;...
    'ModifiedNewton/',true;...
    'QuasiNewton/',true;...
    'ParticleSwarm/',true;...
    'PenaltyBarrierFunctions',true;...
    'QuadraticPrograming/',true;...
    'SequentialQuadraticPrograming/',true;...
    'SteepestDescent/',true;...
    };

for i=1:size(addDirectory,1)
    currentPath = [toolboxFolderPath,addDirectory{i,1}];

    if(addDirectory{i,2})
        % also add subfolders
        totalCurrentPath = genpath(currentPath);
    else
        % add only main folder
        totalCurrentPath = currentPath;
    end

    if(exist(currentPath, 'dir'))
       addpath(totalCurrentPath);
    end
end
clear addDirectory i currentPath totalCurrentPath


%% clear remaining variables
clear toolboxFolderPath  


%% Add SSO-Toolbox
run('./sso-toolbox/setup_sso_toolbox.m');

