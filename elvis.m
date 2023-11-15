function elvisutil = elvis(Z, thetain, gamin, rep, Moment, ...
	guess_un, jump_un, dimf, rho)
thetao = thetain;
   gmm_ = @(thetagam) gmm(Z, rep, Moment, thetao, ...
            dimf, guess_un, jump_un, rho,thetagam);

x0 = [thetain, gamin];
lb = x0-1;
ub = x0+1;
options = optimoptions('fmincon','UseParallel',true);
opt_sol = fmincon(gmm_,x0,[],[],[],[],lb, ...
    ub,[],options);
 thetasoln = opt_sol;
 elvisutil  = thetasoln;
end
