function [Nmodel] = singleexp(p, N0, t)

% This function is given a growth rate parameter p (1) and an initial cell 
% number N0 and returns the corresponding model cell numbers


% model
% N(t) = N_0 * exp(g*t)


Nmodel = N0.* exp((p(1)*t));


end