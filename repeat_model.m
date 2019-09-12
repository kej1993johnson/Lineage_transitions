function [pno_repeat, Enunique, Enreps]= repeat_model(N0,q)

probi(1)=1;
for i = 2:N0
    probi(i) = ((q-i+1)./q);
end

pno_repeat = prod(probi);


Enreps = N0*(1-((1-1/q).^(N0-1)));

Enunique = N0*((1-(1/q)).^(N0-1));




end