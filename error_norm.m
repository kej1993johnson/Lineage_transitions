function error = error_norm(param, sigma, N0, Nf, tf)

P = num2cell(param); 
[g_mu] = deal(P{:}); % our parameters a & b

cell_growth_rates = normrnd(g_mu,sigma,N0,1);
N_cells = ones(N0,1);

for j = 1:length(cell_growth_rates)
    N_cells(j,1)= exp(cell_growth_rates(j).*tf);
end

sum_Ncells= sum(N_cells);

error = (sum_Ncells-Nf).^2;

end