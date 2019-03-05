function [ err ] = fit_mu_var( p, mu_t, var_t,tvec, N0, V0 )
 % function finds error between model and measured mean and variance
 mu_model = mu_fxn(p,tvec, N0);
V_model= V_fxn(p,tvec,N0, V0);
err_mu = mu_model- mu_t;
err_var = V_model- var_t;
err = vertcat(err_mu, err_var);


end

