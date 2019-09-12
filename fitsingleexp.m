function [ err ] = fitsingleexp( param, Nmeas, N0, t)
%This function returns a cost function to be minimized
    g = param(1);
    N_model = singleexp(g, N0,t);
    err = N_model(2:end) - Nmeas;


end