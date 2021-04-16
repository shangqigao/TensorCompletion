function cp = ktensor_approx(x, R, varargin)
%KTENSOR_APPROX Approximation of htensor by ktensor.
%
%   CP = KTENSOR_APPROX(X, R) returns a Tensor Toolbox ktensor of rank R, 
%   which approximates the htensor X. This function requires the Tensor 
%   Toolbox to be installed.
%
%   CP = KTENSOR_APPROX(X, R, varargin) allows to specify more arguments.
%   See CP_ALS.
%
%   This is a wrapper function for the function CP_ALS of the Tensor
%   Toolbox. It was developed and tested with Tensor Toolbox Version 2.4.
%
%   See also CP_ALS, MTTKRP.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

cp = cp_als(x, R, varargin{:});
