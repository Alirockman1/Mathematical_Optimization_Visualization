function fitnessValue = genetic_algorithm_fitness_inverse(objectiveValue,varargin)
%GENETIC_ALGORITHM_FITNESS_INVERSE computes fitness values using inverse 
%   transformation for genetic algorithm optimization. This function converts
%   objective values to fitness values by taking the reciprocal, making
%   smaller objective values correspond to higher fitness values.
%
%	FITNESSVALUE = GENETIC_ALGORITHM_FITNESS_INVERSE(OBJECTIVEVALUE) computes
%	fitness values by taking the inverse of the objective values. This is
%	suitable for minimization problems where smaller objective values should
%	have higher fitness.
%
%	FITNESSVALUE = GENETIC_ALGORITHM_FITNESS_INVERSE(OBJECTIVEVALUE,VARARGIN)
%	allows for additional options to be passed (currently unused but maintained
%	for consistency with other fitness functions).
%
%   Note: This function assumes all objective values are positive. For
%   objective values that can be negative or zero, consider using other
%   fitness transformation methods.
%
%   Input:
%		- OBJECTIVEVALUE : (nPopulation,1) double - objective function values
%		- VARARGIN : additional options (currently unused)
%
%   Output:
%		- FITNESSVALUE : (nPopulation,1) double - fitness values (inverse of objective values)
%
%   Example:
%       objectiveValues = [2; 4; 1; 8];
%       fitnessValues = genetic_algorithm_fitness_inverse(objectiveValues);
%       % Result: fitnessValues = [0.5; 0.25; 1.0; 0.125]
%
%   See also GENETIC_ALGORITHM_SELECTION_ROULETTE_WHEEL, OPTIMIZATION_GENETIC_ALGORITHM
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
	fitnessValue = 1./objectiveValue;
end