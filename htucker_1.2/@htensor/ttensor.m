function x_tucker = ttensor(x)
%TTENSOR Convert htensor to a Tensor Toolbox ttensor.
%
%   Y = TTENSOR(X) returns a Tucker tensor in the Tensor Toolbox
%   format. This function requires the Tensor Toolbox to be installed.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

% Setting U{ii} to the identity matrix reduces x to its core
% tensor.
U = x.U(x.dim2ind);

rk = rank(x);
for ii=find(x.is_leaf)
  x.U{ii} = eye(rk(ii));
end

% Calculate full core tensor
core = tensor(full(x));

% Construct ttensor
x_tucker = ttensor(core, U);