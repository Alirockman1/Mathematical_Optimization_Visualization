function [inverseQuasiHessian, inverseQuasiHessianData] = quasi_newton_bfgs(currentDesign, currentGradient, nextDesign, nextGradient, currentInverseQuasiHessian, varargin)
%QUASI_NEWTON_BFGS Computes the Hessian matrix using the BFGS update and its inverse using the Sherman-Morrison formula.
%   This function approximates the Hessian using the BFGS update and maintains
%   the inverse Hessian using the Sherman-Morrison formula.
%
%   [HESSIAN, INVERSEHESSIAN] = QUASI_NEWTON_BFGS(OBJECTIVEFUNCTION, DESIGNSAMPLE)
%   computes the initial Hessian approximation and updates it using BFGS.
%
%   Inputs:
%       - OBJECTIVEFUNCTION : function handle
%       - DESIGNSAMPLE : (nDimension, 1) double (initial design point)
%
%   Outputs:
%       - HESSIAN : (nDimension, nDimension) double (updated Hessian)
%       - INVERSEHESSIAN : (nDimension, nDimension) double (inverse Hessian)
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

    if(isempty(currentGradient) || isempty(nextGradient))
        inverseQuasiHessian = currentInverseQuasiHessian;
        inverseQuasiHessianData = [];
        return;
    end

    % Compute Pertubations
    changeOfGradient = nextGradient - currentGradient;
    changeOfDesign = nextDesign - currentDesign;

    % transform to column vectors for compatibility with usual mathematical notation
    changeOfGradient = changeOfGradient(:);
    changeOfDesign = changeOfDesign(:);
    
    % Terms to compile the Sherman Morrison formula
    firstTerm  = currentInverseQuasiHessian;
    secondTerm = (1 + ((changeOfGradient'*currentInverseQuasiHessian*changeOfGradient)/(changeOfDesign'*changeOfGradient)))*((changeOfDesign*changeOfDesign')/(changeOfDesign'*changeOfGradient));
    thirdTerm  = ((changeOfDesign*changeOfGradient'*currentInverseQuasiHessian + currentInverseQuasiHessian*changeOfGradient*changeOfDesign')/(changeOfDesign'*changeOfGradient));

    % Sherman-Morrison formula
    inverseQuasiHessian = firstTerm + secondTerm - thirdTerm;

    % Prepare output data if requested
    isOutputHessianData = (nargout >= 2);
    if(isOutputHessianData)
        inverseQuasiHessianData = struct(...
            'CurrentInverseQuasiHessian', currentInverseQuasiHessian, ...
            'NextInverseQuasiHessian', inverseQuasiHessian, ...
            'DeltaGradient', changeOfGradient, ...
            'DeltaDesign', changeOfDesign, ...
            'BfgsFirstTerm', firstTerm, ...
            'BfgsSecondTerm', secondTerm, ...
            'BfgsThirdTerm', thirdTerm);
    end
end
