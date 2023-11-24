%% Metropolis Hastings Algorithm
function avg_mom = mom_ave(Z, rep, Moment, thetain, dimf, guess_un, ...
    jump_un, rho, thetagam)  
 		theta = thetagam(1:length(thetain));
        n = length(Z);
 		gam = thetagam((length(thetain)+1):(length(thetain)+dimf));
 		r = -rep(1)+1;
 		%Store average as a for accepted draws
 		avg = zeros(dimf, n);
 		%Store accepted values of the vector of unobservables
 		unmat = zeros(n, rep(2));
 		%Initial Guess for Unobservables
 		un = guess_un;
    for i = r:rep(2)
			%Proposal distribution 
			try_un = jump_un;
			%Acceptance Probability
			try_dens = exp((gam*(Moment(try_un, Z, theta)))./ ...
                (gam*(Moment(un, Z, theta)))).* ...
			(rho./rho);
			%Create a set of indices for acceptance rule
			aindex = find(rand([1,n]) < try_dens);
			%Update unobservables that satisfy acceptance rule
			un(aindex) = try_un(aindex);
        if i>0 
				unmat(:,i) = un;
				avg = (avg+Moment(un, Z, theta))/rep(2);
        end
    end
        avg(isnan(avg))=0;
 		avg_mom = avg;
    end