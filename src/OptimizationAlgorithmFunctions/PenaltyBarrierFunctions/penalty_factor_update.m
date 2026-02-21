function penaltyFactor = penalty_factor_update(initialPenalty,currentPenalty,iteration)
%PENALTY_FACTOR_UPDATE updates the penalty factor for penalty and barrier methods
%   in constrained optimization problems.
%
%	PENALTYFACTOR = PENALTY_FACTOR_UPDATE(INITIALPENALTY,CURRENTPENALTY,ITERATION) 
%   updates the penalty factor by reducing the current penalty factor by half.
%   This is a simple reduction strategy commonly used in penalty methods where
%   the penalty factor is progressively reduced to improve convergence behavior.
%
%	The function implements a basic penalty factor update rule where:
%   - The penalty factor is halved at each iteration
%   - This reduction helps balance constraint satisfaction with objective optimization
%   - The initial penalty and iteration parameters are available for more complex
%     update strategies but are not used in this basic implementation
%
%   Input:
%		- INITIALPENALTY : double - initial penalty factor value
%		- CURRENTPENALTY : double - current penalty factor value  
%		- ITERATION : integer - current iteration number
%
%   Output:
%		- PENALTYFACTOR : double - updated penalty factor (currentPenalty/2)
%
%   Example:
%       % Basic usage in penalty method
%       initialPenalty = 10;
%       currentPenalty = 5;
%       iteration = 1;
%       newPenalty = penalty_factor_update(initialPenalty, currentPenalty, iteration);
%       % newPenalty will be 2.5
%
%   See also OPTIMIZATION_BARRIER_PENALTY, PENALTY_EQUALITY_FUNCTION_SQUARE,
%   PENALTY_INEQUALITY_FUNCTION_SQUARE, BARRIER_FUNCTION_INVERSE,
%   BARRIER_FUNCTION_NATURAL_LOG
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

	penaltyFactor = currentPenalty/2;
end