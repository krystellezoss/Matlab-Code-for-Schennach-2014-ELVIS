%% Routine implementing ELVIS in Matlab

% general naming conventions
% un:unobservables
% theta: parameter vector
% nuis: nuisance parameters (profiled over)
% gam: tilting parameter gamma (see paper for details)

% Calculate average of moment function over the unobservables for each individual
% Moment: pointer to moment function
% theta: parameter vector, including nuisance parameters
% gam: tilting parameter gamma
% rep: array controling replications:
% 	rep[1]: number of equilibration steps
% 	rep[2]: number of averaging steps
% Z: Observables
% dimf: number of moment equations. Needs to match the dimension of Moment
%thetain: Initial values of theta to search over
%gamin: Initial values of tilting parameter, gamma to search over (must match the dimension of dimf)


%% ELVIS
%moments = @(U, Z, theta) ...
%((Z(1,:)+U-Z(:,2)*theta(1)).*Z(2,:));

%Sample Size
n = 250;
res = 1;
%True value of Theta[1]
thetao = [1, 0.2];
x = normrnd(0, 1, [1,n]);
	dy = 0.5*normrnd(0, 1, [1,n]);
	yt = thetao(1)*x+dy;
	y = res*floor(yt/res);
	Z = [y; x];
    U = rand([1,n]);


guess_un = rand([1,n]);

jump_un = rand([1,n]);

rho= ones(1, n);

%% MCMC Rules 
%Burn in and integration replications
rep = [50,500];
%Number of MC replications
nreps = 50;


%% DGP
%Interval Valued Data Example
%Moment Function
 Moment = @(U, Z, theta) ...
	([Z(1,:)+U-Z(2,:)*theta(1); (Z(1,:)+U-Z(2,:)*theta(1)).*Z(2,:); ... 
	((Z(1,:)+U-Z(2,:)*theta(1)).^2)-theta(2); (Z(2,:).^2).* ...
    (((Z(1,:)+U-Z(2,:)*theta(1)).^2)-theta(2))]);


%Matrix to store objective function values over grid
soln = zeros(nreps, 6);
tic
for j = 1:nreps	
    size_ = size(Moment(U, Z, thetao));
	dimf = size_(1);
    gamo = ones(1,dimf);
	soln(j,:) = elvis(Z, thetao, gamo, rep, Moment, ...
	guess_un, jump_un, dimf, rho);
end
mean_theta = mean(soln(:, 1:2), 1);
se_theta = std(soln(:, 1:2));
toc
