%% Define MM Objective Function


function GMM_objective = gmm(Z, rep, Moment, thetao, ...
            dimf, guess_un, jump_un, rho,thetagam)

 		mom = mom_ave(Z, rep, Moment, thetao, ...
            dimf, guess_un, jump_un, rho, thetagam); 

 		obj = length(mom)*mean(mom, 1)*(cov(mom)\ ...
        transpose(mean(mom,1)));
 		GMM_objective = -obj;
end
