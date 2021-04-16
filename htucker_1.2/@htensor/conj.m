function x = conj(x)
%CONJ Complex conjugate of htensor.
%
%  CONJ(x) returns the complex conjugate of htensor x.
%
%  See also HTENSOR/PERMUTE

% htucker toolbox
% C. Tobler and D. Kressner, EPF Lausanne
% FreeBSD License, see COPYRIGHT.txt

x_is_leaf = x.is_leaf;

for ii=1:x.nr_nodes
  if(x_is_leaf(ii))
    x.U{ii} = conj(x.U{ii});
  else
    x.B{ii} = conj(x.B{ii});
  end
end