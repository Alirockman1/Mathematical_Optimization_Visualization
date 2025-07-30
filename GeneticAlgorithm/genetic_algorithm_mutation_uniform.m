function children = genetic_algorithm_mutation_uniform(parents,designSpace,varargin)
%GENETIC_ALGORITHM_MUTATION_UNIFORM performs uniform mutation operation
%   for genetic algorithm optimization. This function mutates selected
%   individuals by replacing one randomly chosen design variable with a
%   uniformly distributed random value within the design space bounds.
%
%	CHILDREN = GENETIC_ALGORITHM_MUTATION_UNIFORM(PARENTS,DESIGNSPACE) performs
%	uniform mutation on the PARENTS population within the specified DESIGNSPACE.
%	For each parent, one design variable is randomly selected and replaced with
%	a uniformly distributed random value between the lower and upper bounds.
%
%	CHILDREN = GENETIC_ALGORITHM_MUTATION_UNIFORM(PARENTS,DESIGNSPACE,VARARGIN)
%	allows for additional options to be passed (currently unused but maintained
%	for consistency with other mutation functions).
%
%   The mutation operation:
%   1. Randomly selects one design variable for each parent
%   2. Generates a uniform random value within the design space bounds
%   3. Replaces the selected variable with the random value
%
%   Input:
%		- PARENTS : (nParents,nDesignVariable) double - parent population to mutate
%		- DESIGNSPACE : (2,nDesignVariable) double - design space bounds [lower; upper]
%		- VARARGIN : additional options (currently unused)
%
%   Output:
%		- CHILDREN : (nParents,nDesignVariable) double - mutated population
%
%   Example:
%       parents = [1.5, 2.3; 0.8, 1.9];
%       designSpace = [0, 0; 3, 4];
%       children = genetic_algorithm_mutation_uniform(parents, designSpace);
%       % One variable in each parent will be randomly replaced
%
%   See also GENETIC_ALGORITHM_MUTATION_BOUNDARY, OPTIMIZATION_GENETIC_ALGORITHM
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

	populationSize = size(parents,1);
	nDimension = size(parents,2);

	boundaryChange = randsample(nDimension,populationSize,true);
	uniformDistribution = rand(populationSize,1);

	children = parents;
	children(:,boundaryChange) = designSpace(1,boundaryChange) + uniformDistribution.*(designSpace(2,boundaryChange)-designSpace(1,boundaryChange));
end