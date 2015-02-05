function [a,b] = solve_cbp(A,B,y,lambda, options)

% solve_cbp - solve the C-BP using CVX
%
%   [a,b] = solve_cbp(A,B,y,lambda);
%
% for lambda=0, solve
%       min |a|_1  s.t. A*a+B*b=y   and  a>=0,  |b|<=a
% for lambda>0, solve
%       min |a|_1  + 1/(2*lambda) * |A*a+B*b - y|^2   and  a>=0, |b|<=a
%
% In particular, for B=0 (or B=[]), solves positive BP.
%
%   You must have CVX installed.
%
%   Copyright (c) 2014 Gabriel Peyre

N = size(A,2);

if isempty(B)
    % Positive LASSO
    B = A*0;
end

options.null = 0;
verb = getoptions(options, 'verb', 'quiet');
precis = getoptions(options, 'precision', 'high'); % or 'best'
    
cvx_solver SeDuMi % sdpt3 %  %
cvx_begin(verb);
cvx_precision(precis);
variable a(N);
variable b(N);
abs(b)<=a;  % a is automatically positive
if lambda==0
    A*a + B*b==y;
    minimize( sum(a) )
else
    minimize( 1/2*sum(pow_abs(A*a + B*b-y,2)) + lambda*sum(a) )
end
cvx_end

end
