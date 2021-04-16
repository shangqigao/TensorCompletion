function example_cancellation()
%EXAMPLE_CANCELLATION Cancellation in tan(x) + 1/x - tan(x).
%
%   Calculate elementwise function tan(x) + 1/x - tan(x), where x is a 3rd
%   order tensor containing the sum
%         xi_1+xi_2+xi_3  for  xi_i = linspace(0, 1, 101).
%   Let T denote truncation to HTD. Then we compare the following two
%   methods:
%      Successive truncation: T( T( tan(x) + 1/x ) - tan(x) ),
%      Direct computation:    T( tan(x) + 1/x - tan(x) ).
%
%   See also EXAMPLE_CANCELLATION2.

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

close all;

fprintf(['Calculate tan(x) + 1/x - tan(x), using successive truncation\n' ...
         'and directly.\n']);

d = 3;
n = 100;
opts.max_rank = 20;

t = linspace(0, 1, n+1)';
t = t(2:end);

% Construct x1+x2+x3
sum_t = cell(1, d);
for ii=1:d
  sum_t{ii} = ones(n, d);
  sum_t{ii}(:, ii) = t;
end
sum_full = full(htensor(sum_t));

% Construct 1/(x1+x2+x3)
inv_sum_full = 1./sum_full;
inv_sum_htd = truncate(inv_sum_full, opts);

% Construct tan(x1+x2+x3)
tan_full = tan(sum_full);
tan_htd = truncate(tan_full, opts);


sum_exact = tan_htd + inv_sum_htd - tan_htd;

sum_trunc_standard = truncate(tan_htd + inv_sum_htd - tan_htd, ...
			      opts);

sum_trunc_sum = htensor.truncate_sum({tan_htd, inv_sum_htd, -tan_htd}, ...
			      opts);

sum_succ_trunc = truncate(tan_htd + inv_sum_htd, opts);
sum_succ_trunc = truncate(sum_succ_trunc - tan_htd, opts);

err_std_trunc  = norm(sum_exact - sum_trunc_standard)/norm(sum_exact)
err_trunc_sum  = norm(sum_exact - sum_trunc_sum)/norm(sum_exact)
err_succ_trunc = norm(sum_exact - sum_succ_trunc)/norm(sum_exact)

fprintf(['Note that in this (constructed) example, the truncation error\n' ...
         'leads to massive cancellation in the case of successive\n'...
         'truncation.\n']);