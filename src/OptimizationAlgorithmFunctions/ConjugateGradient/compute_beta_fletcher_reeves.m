function beta = compute_beta_fletcher_reeves(gradientCurrent, gradientPrevious)
%COMPUTE_BETA_FLETCHER_REEVES Fletcher-Reeves beta parameter for conjugate gradient.
%
%   BETA = COMPUTE_BETA_FLETCHER_REEVES(GRADCURRENT, GRADPREVIOUS) computes
%   the Fletcher-Reeves conjugate gradient parameter β according to:
%       β = (∇f_k+1ᵀ∇f_k+1) / (∇f_kᵀ∇f_k)
%
%   INPUTS:
%       GRADCURRENT   - Gradient vector at current iteration (k+1)
%       GRADPREVIOUS  - Gradient vector at previous iteration (k)
%
%   OUTPUTS:
%       BETA          - Fletcher-Reeves beta parameter
%
%   FUNCTIONALITY:
%       - Computes the Fletcher-Reeves formula for conjugate gradient methods
%       - Automatically handles division by zero cases
%       - Uses efficient norm computations
%
%   NOTES:
%       - Returns β = 0 when denominator would be zero
%       - Uses Euclidean norm (L2-norm) for gradient magnitudes
%       - Numerically stable implementation
%
%   REFERENCES:
%       [1] Fletcher, R., Reeves, C.M. (1964). "Function minimization by
%           conjugate gradients". The Computer Journal 7 (2): 149-154.
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

    beta = norm(gradientCurrent)^2 / norm(gradientPrevious)^2; 
    
    % Avoid division by zero
    if(isnan(beta) || isinf(beta))
        beta = 0;
    end
end