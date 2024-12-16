function LLH=Ranking_LLH(ypred,yobs)
    %% LH & LLH
    for i = 1:size(ypred,2)
       Yhat = ypred(:,i);
        Res = log(yobs./981) - log(Yhat./981);
        sigma = std(Res);
    
    
        xi = log(yobs);
        mu = log(Yhat);
        sigma(1,i) = sqrt(var(Res));
    
        z = (xi-mu)/sigma(1,i);
    %%
        z1 = abs(z);
        LHi = 2*1/2*(erf(inf)-erf(z1/sqrt(2)));
        LH_index(1,i) = median(LHi); format short;
        meanZ(1,i)= mean(z);
    
        LLHi= -log2(normpdf(xi,(mu),sigma(1,i)));
        LLH(1,i) = mean(LLHi);
    %%
end