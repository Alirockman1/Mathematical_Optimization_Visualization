function children = genetic_algorithm_crossover_arithmetic(parents,fitnessValue,varargin)
%GENETIC_ALGORITHM_CROSSOVER_ARITHMETIC performs arithmetic crossover operation
%   for genetic algorithm optimization. This function creates offspring by
%   computing weighted averages of two parents, where the weights are
%   randomly generated and complementary (sum to 1).
%
%	CHILDREN = GENETIC_ALGORITHM_CROSSOVER_ARITHMETIC(PARENTS,FITNESSVALUE) performs
%	arithmetic crossover on two PARENTS. The FITNESS values are passed for
%	consistency with other crossover functions but are not used in the
%	arithmetic crossover operation.
%
%	CHILDREN = GENETIC_ALGORITHM_CROSSOVER_ARITHMETIC(PARENTS,FITNESSVALUE,VARARGIN)
%	allows for additional options to be passed (currently unused but maintained
%	for consistency with other crossover functions).
%
%   The crossover operation:
%   1. Generates a random weight factor α ∈ [0,1]
%   2. First child = α * parent1 + (1-α) * parent2
%   3. Second child = (1-α) * parent1 + α * parent2
%
%   This method ensures that both children lie within the convex hull of
%   the two parents, promoting intermediate solutions.
%
%   Input:
%		- PARENTS : (2,nDesignVariable) double - two parent individuals
%		- FITNESSVALUE : (2,1) double - fitness values (not used in arithmetic crossover)
%		- VARARGIN : additional options (currently unused)
%
%   Output:
%		- CHILDREN : (2,nDesignVariable) double - two offspring individuals
%
%   Example:
%       parents = [1.0, 2.0; 3.0, 4.0];
%       fitnessValue = [0.5; 0.8];  % Not used in arithmetic crossover
%       children = genetic_algorithm_crossover_arithmetic(parents, fitnessValue);
%       % If randomNumber = 0.3:
%       % children(1,:) = 0.3*[1,2] + 0.7*[3,4] = [2.4, 3.4]
%       % children(2,:) = 0.7*[1,2] + 0.3*[3,4] = [1.6, 2.6]
%
%   See also GENETIC_ALGORITHM_CROSSOVER_HEURISTIC, GENETIC_ALGORITHM_CROSSOVER_VARIABLE_EXCHANGE, OPTIMIZATION_GENETIC_ALGORITHM
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

	randomNumber = rand;
	children(1,:) = randomNumber.*parents(1,:) + (1-randomNumber).*parents(2,:);
	children(2,:) = (1-randomNumber).*parents(1,:) + randomNumber.*parents(2,:);
end