function selectionIndex = genetic_algorithm_selection_roulette_wheel(fitnessValue,varargin)
%GENETIC_ALGORITHM_SELECTION_ROULETTE_WHEEL performs roulette wheel selection
%   for genetic algorithm optimization. This function selects individuals
%   from a population based on their fitness values, where individuals with
%   higher fitness have a higher probability of being selected.
%
%	SELECTIONINDEX = GENETIC_ALGORITHM_SELECTION_ROULETTE_WHEEL(FITNESSVALUE)
%	performs roulette wheel selection based on the FITNESS values of the
%	population. Each individual's selection probability is proportional to
%	its fitness value relative to the total population fitness.
%
%	SELECTIONINDEX = GENETIC_ALGORITHM_SELECTION_ROULETTE_WHEEL(FITNESSVALUE,VARARGIN)
%	allows for additional options to be passed (currently unused but maintained
%	for consistency with other selection functions).
%
%   The selection operation:
%   1. Calculates selection probability for each individual as:
%      P(i) = fitness(i) / sum(fitness)
%   2. Performs weighted random sampling with replacement
%   3. Returns indices of selected individuals
%
%   This method ensures that fitter individuals have higher chances of
%   being selected while still allowing less fit individuals to potentially
%   contribute to the next generation.
%
%   Input:
%		- FITNESSVALUE : (nPopulation,1) double - fitness values for population
%		- VARARGIN : additional options (currently unused)
%
%   Output:
%		- SELECTIONINDEX : (nPopulation,1) double - indices of selected individuals
%
%   Example:
%       fitnessValue = [0.1; 0.3; 0.6];  % Population of 3 individuals
%       selectionIndex = genetic_algorithm_selection_roulette_wheel(fitnessValue);
%       % Individual 3 has 60% chance, Individual 2 has 30% chance, Individual 1 has 10% chance
%       % selectionIndex might be [3; 2; 3] (with replacement)
%
%   Note: This function assumes all fitness values are non-negative. For
%   negative fitness values, consider fitness scaling or other selection methods.
%
%   See also GENETIC_ALGORITHM_FITNESS_INVERSE, OPTIMIZATION_GENETIC_ALGORITHM
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
	probabilityOfSelection = fitnessValue / sum(fitnessValue);
	populationSize = size(fitnessValue,1);
	
    selectionIndex = randsample(populationSize, populationSize, true, probabilityOfSelection);
end