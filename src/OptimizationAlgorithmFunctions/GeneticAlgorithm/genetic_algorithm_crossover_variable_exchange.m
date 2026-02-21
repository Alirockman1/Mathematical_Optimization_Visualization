function children = genetic_algorithm_crossover_variable_exchange(parents,fitnessValue,varargin)
%GENETIC_ALGORITHM_CROSSOVER_VARIABLE_EXCHANGE performs variable exchange 
%   crossover operation for genetic algorithm optimization. This function 
%   creates offspring by exchanging design variables between two parents 
%   based on a randomly generated binary mask.
%
%	CHILDREN = GENETIC_ALGORITHM_CROSSOVER_VARIABLE_EXCHANGE(PARENTS,FITNESSVALUE) 
%	performs variable exchange crossover on two PARENTS. The FITNESS values 
%	are passed for consistency with other crossover functions but are not 
%	used in the variable exchange operation.
%
%	CHILDREN = GENETIC_ALGORITHM_CROSSOVER_VARIABLE_EXCHANGE(PARENTS,FITNESSVALUE,VARARGIN)
%	allows for additional options to be passed (currently unused but maintained
%	for consistency with other crossover functions).
%
%   The crossover operation:
%   1. Generates a random binary mask for each design variable
%   2. For each variable where mask = true:
%      - First child gets the variable from second parent
%      - Second child gets the variable from first parent
%   3. For each variable where mask = false:
%      - Children retain variables from their respective parents
%
%   This method allows for discrete exchange of design variables while
%   maintaining the structure of both parents.
%
%   Input:
%		- PARENTS : (2,nDesignVariable) double - two parent individuals
%		- FITNESSVALUE : (2,1) double - fitness values (not used in variable exchange)
%		- VARARGIN : additional options (currently unused)
%
%   Output:
%		- CHILDREN : (2,nDesignVariable) double - two offspring individuals
%
%   Example:
%       parents = [1.0, 2.0, 3.0; 4.0, 5.0, 6.0];
%       fitnessValue = [0.5; 0.8];  % Not used in variable exchange
%       children = genetic_algorithm_crossover_variable_exchange(parents, fitnessValue);
%       % If mask = [true, false, true]:
%       % children(1,:) = [4.0, 2.0, 6.0] (exchanged variables 1 and 3)
%       % children(2,:) = [1.0, 5.0, 3.0] (exchanged variables 1 and 3)
%
%   See also GENETIC_ALGORITHM_CROSSOVER_ARITHMETIC, GENETIC_ALGORITHM_CROSSOVER_HEURISTIC, OPTIMIZATION_GENETIC_ALGORITHM
%   
%   Copyright 2025 Eduardo Rodrigues Della Noce
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
	children = parents;
	nDimension = size(parents,2);

	% Generate random crossover mask
    mask = (rand(1,nDimension)>0.5);
    children(1,mask) = parents(2, mask);
    children(2,mask) = parents(1, mask);
end