function [x0,Slist] = gen_sparse_signals(slist,n,q)

% gen_sparse_signals - genetate sparse signal
%
%   [x0,Slist] = gen_sparse_signals(slist,n,q);
%
% Generate a list of q signals within a range of sparsity   
%   x0(:,i) is Slist(i) sparse.
%
%   Copyright (c) 2015 Gabriel Peyre

% List of sparsity of each signal
Slist = slist(mod(0:q-1,length(slist))+1);

% Genetate signals so that |x0(:,j)| has sparsity |Slist(i)|.
U = rand(n,q);
v = sort(U);
v = v( (0:q-1)*n + Slist );
x0 = double( U<=repmat( v, [n 1] ) );
x0 = x0.*sign(randn(size(x0)));