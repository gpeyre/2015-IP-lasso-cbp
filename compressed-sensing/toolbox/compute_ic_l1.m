function ic = compute_ic_l1(X,Phi)

% compute_ic_l1 - compute Fuchs criterion
%
%   ic = compute_ic_l1(X,Phi);
%
% Copyright (c) Gabriel Peyre 2014

epsilon = 1e-8;
I = find(abs(X)>epsilon);
J = find(abs(X)<=epsilon);

PhiI = Phi(:,I);
PhiJ = Phi(:,J);

ic = norm( PhiJ' * pinv(PhiI)' * sign(X(I)) , 'inf' );
