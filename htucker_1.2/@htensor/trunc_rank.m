function [k, err, success] = trunc_rank(s, opts)
%TRUNC_RANK Return rank according to user-specified parameters.
%
%   K = TRUNC_RANK(S, OPTS) returns the rank K to truncate the singular
%   values in the array S according to the requirements in OPTS (see
%   TRUNCATE). This is an internal function used in all truncation
%   functions.
%
%   [K, ERR, SUCCESS] = TRUNC_RANK(S, OPTS) also returns the error ERR
%   and a flag SUCCESS signalling whether all requirements specified in
%   OPTS are satisfied.
%
%   See also: HTENSOR, TRUNCATE.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

if(nargin ~= 2)
  error('Requires 2 arguments.')
end

if(~isnumeric(s) || ~isvector(s))
  error('First argument must be a vector.');
end

if(~isa(opts, 'struct') || ~isfield(opts, 'max_rank') )
  error(['Second argument must be a MATLAB struct with ' ...
	 ' field max_rank (and optionally other fields).']);
end

% When truncating at k, error in Frobenius norm is s_sum(k+1) with
% s_sum(k+1) = norm(s(k+1:end))
s_sum = sqrt(cumsum(s(end:-1:1).^2));
s_sum = s_sum(end:-1:1);

% Calculate necessary rank to satisfy relative eps
if(~isfield(opts, 'rel_eps'))
  k_rel = [];
else
  k_rel = find(s_sum < opts.rel_eps*norm(s), 1, 'first');
  k_rel = k_rel-1;
  if(numel(k_rel) == 0)
    k_rel = length(s);
  end
end

% Calculate necessary rank to satisfy absolute eps
if(~isfield(opts, 'abs_eps'))
  k_abs = [];
else
  k_abs = find(s_sum < opts.abs_eps, 1, 'first');
  k_abs = k_abs-1;
  if(numel(k_abs) == 0)
    k_abs = length(s);
  end
end

opts.max_rank = min([opts.max_rank, numel(s)]);

% Calculate necessary rank to satisfy absolute and relative eps
k_eps = max([k_rel, k_abs]);

% Use k_eps if <= maximal rank, otherwise use maximal rank
k = min([k_eps, opts.max_rank]);
k = max(k, 1);

% Read error from s_sum
s_sum = [s_sum; 0];
err = s_sum(k+1);

% Warning flag when relative or absolute eps are not satisfied.
if(k < k_eps)
  success = false;
else
  success = true;
end

% For debugging / illustration: Plot the values s, rel_eps, abs_eps
% and the truncation rank k.
if(isfield(opts, 'plot_sv') && opts.plot_sv)
  semilogy(0:length(s_sum)-1, s_sum, 'b-x');
  xlim([0 length(s)])
  hold on;
  if(isfield(opts, 'rel_eps'))
    semilogy([0 length(s)], opts.rel_eps*norm(s)*[1 1], 'r');
  end
  if(isfield(opts, 'abs_eps'))
    semilogy([0 length(s)], opts.abs_eps*[1 1], 'g');
  end
  if(opts.max_rank <= length(s))
    semilogy(opts.max_rank*[1 1], ylim(), 'k');
  end
  semilogy(k, s_sum(k+1), 'ro');
  
  legend('error decay', 'rel. eps', 'abs. eps', 'max. rank') 
  drawnow;
end
