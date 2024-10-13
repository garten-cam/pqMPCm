% Script for testing the pqMPCm algorithm with the simple bioreactor.
% pqMPCm : Model Predictive Control of a system based on a pqEDMD
% decompostion for Matlab.
% 
% The script assumes that the pqEDMDm package is also in the Matlab path.
%
%-------------------------------------------------------------------------%
%------------ The Experiments --------------------------------------------%
%-------------------------------------------------------------------------%
% We will assume a "realistic" scenario where the data comes from
% experiments of the process where the initial conditions of substrate
% concentration and bionmass are small, and in operating conditions where
% the biomass thrives, this excludes "failed" experiments.
%
clear variables
rng(1) % for consistency of simulation
%
% Define the parameters for the simulation
% num_ics = 10;
% From a unique initial condition
ics = [(0.15-0.1)*rand(1)+0.1, (0.2-0.1)*rand(1)+0.1];
%
tfin = 950; % minus one the real tfin is tfin-1po
%
sigma = 0.1;
sampling_factor = .2;
n_points = floor(sampling_factor*tfin)+1; 
t = linspace(0,tfin,n_points);
% Define a profile for the input
% Number of dillution rates for the allotted time
% n_drs = 4;
% drts = (0.2-0.01)*rand(1,n_drs)+0.01;
% dr_rfs = [0.10, 0.273, 0.010, 0.2731, 0.09];
dr_rfs = [0.15, 0.273, 0.010, 0.273, 0.09, 0.26, 0.09];
% dr_rfs = [0.150, 0.273, 0.010, 0.24, 0.09];
% dr_rfs = [0.02, 0.1, 0.26, 0.05, 0.15];
% divide the time into n_drs segments
ts = floor(linspace(1,n_points,length(dr_rfs)+1)); % Time segments
tins = t(ts(1:end-1)); % Initial times for the segments
tens = t(ts(2:end)); % End or final time for the segments
%
% Define function for the input
step_gen = @(t,tis,tes,frs)sum( ...
	cell2mat(arrayfun(@(ti,te,dr){(t>=ti).*(t<te)*dr},tis,tes,frs)'), ...
	1);
in = @(t)step_gen(t,tins,tens,dr_rfs);
%
%
% Run the simulation
odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
[exp.t, exp.yd] = ode23s(@(t,x)bioODE(t,x, ...
	in),... there is a different mu_max per experiment
	linspace(0,tfin-1,n_points), ...
	ics, ...
	odeSettings);
exp.y = exp.yd + normrnd(0,sigma,size(exp.yd));
exp.u = in(exp.t')';

% Plot the experiment
figure(1)
clf
tiledlayout(2,1,TileSpacing="tight")
nexttile
plot(exp.t,exp.y)
% xlabel('t',Interpreter='latex')
ylabel('$x$',Interpreter='latex')
title('Experimental Data')
legend('Biomass','Substrate')
nexttile
plot(exp.t,exp.u)
ylabel('$D$',Interpreter='latex')
xlabel('t[h]',Interpreter='latex')

%
%-------------------------------------------------------------------------%
%------------ Identification ---------------------------------------------%
%-------------------------------------------------------------------------%
% Divide take the first segments as training and the last for testing
% Normalize the experiments
% exp_nstch = normalize_data(exp_stch,[-1,1]);
exptr = struct('y',exp.y(ts(1):ts(end-1),:),...
	't',exp.t(ts(1):ts(end-1),:),...
	'u',exp.u(ts(1):ts(end-1),:));
expts = struct('y',exp.y(ts(end-1):ts(end),:),...
	't', exp.t(ts(end-1):ts(end),:),...
	'u', exp.u(ts(end-1):ts(end),:),...
	'yd',exp.yd(ts(end-1):ts(end),:));
pen_edmd_stch = pqEDMDm( ...
	p = [2 3 4 5],...
	q = [0.5 1 1.5 2 2.5],...
	observable = @legendreObservable,...
	dyn_dcp = @svdDecomposition...
	);
dcps = pen_edmd_stch.fit(exptr);
%
% Calculate the error
st_err = arrayfun(@(dcp)dcp.error(expts),dcps);
% whos the best
[st_min, st_bt] = min(st_err);
% extract the best
dcp = dcps(st_bt);
% Make the prediction from the test set
pred = dcp.pred_from_test(expts);

%%
% Plot the result
st_fig = figure(2);
% instead of a phase plane, use a grid
clf
tiledlayout(3,1,TileSpacing="tight")
trans = .7;
nexttile([2 1])
hold on
ts_p   = plot(exp.t,exp.yd,'r',LineWidth=.5);
st_trpx = plot(exptr.t,exptr.y(:,1),Color=[0 .4 .7 trans],LineWidth=.5);
st_trps = plot(exptr.t,exptr.y(:,2),Color=[.4 .6 .1 trans],LineWidth=.5);
st_tspx = plot(expts.t,expts.y(:,1),Color=[.3 .7 .9 trans],LineWidth=.5);
st_tsps = plot(expts.t,expts.y(:,2),Color=[.9 .6 .1 trans],LineWidth=.5);
st_pqp = plot(expts.t,pred.y,'-.k',LineWidth=.5);
ylabel('$x$ [g/L]',Interpreter='latex')
legend([st_trpx,st_trps,...
	      st_tspx,st_tsps,...
				ts_p(1),st_pqp(1)],...
				{'$x_1$ tr','$x_2$ tr',...
				 '$x_1$ ts','$x_2$ ts',...
				 'truth', 'appx'}, ...
				 Interpreter="latex",Location="northeast",NumColumns=2)
title({"pqEDMD ",  ...
       " $p$="+num2str(dcp.obs.p) + ...
			 " $q$="+num2str(dcp.obs.q) + ...
			 " $\epsilon$="+num2str(st_min) + ...
			 " $\Delta$t="+num2str(t(2)-t(1)) + "[h]"},Interpreter="latex")
nexttile
plot(exp.t, exp.u,Color=[0 .4 .7]);
xlabel('$t$ [h]',Interpreter='latex')
ylabel('$D$',Interpreter='latex')
legend('Dulution rate',Interpreter="latex",Location="southwest")
set(gcf,'PaperPosition', [0, 0, 16, 8])
% saveas(st_fig,"~/Documents/BioReactorConf/figures/pqEDMD_fast_high.fig")
% saveas(st_fig,"~/Documents/BioReactorConf/figures/pqEDMD_fast_high.eps", "epsc")
%%
%-------------------------------------------------------------------------%
%---------- Control ------------------------------------------------------%
%-------------------------------------------------------------------------%
% The workflow should be:
% 
% 1.	Create the matrices object
% 1.1 dcp: is the decomposition object
% 1.2 Np:  is the Prediction horizon
Np = 7;
Nc = 5;
% 1.3 Qy:  is the output/state weight
% 1.4 Qu:	 is the input weight
% 1.5 S:   is the terminal weight 
mpc_mat = pqMPCu_mat(dcp, ... Decompostition objedct
	Np, ... Np
	Nc, ...
	[1 0; 0 0], ... Qy only interested in biomass
	.2*eye(dcp.m), ... Qu 
	1*[1 0; 0 0], ... S The same, only biomass
	"Obs"); % 'Cob' for C matrix of the observables, or 'Obs' for the actual observables

% Main loop for the MPC...
ctr_fin = 200; % final time
% Number of points according to the sampling factor
ctr_points = floor(sampling_factor*ctr_fin)+1;
% time vector for the control
ctr_t = linspace(0,ctr_fin,ctr_points);
Dt = ctr_fin/(ctr_points-1);
% Setpoints for the Biomass
bm_frs = [0.8 1.2 0.7 1.2];
%
ctr_ts = floor(linspace(1,ctr_points, length(bm_frs)+1));
ctr_tins = ctr_t(ctr_ts(1:end-1));
ctr_tens = ctr_t(ctr_ts(2:end));
% reference
bio_ref = @(t)step_gen(t,ctr_tins,ctr_tens,bm_frs);

% Preallocate the output choose a random initial condition
ctr_sim = struct('y',zeros(ctr_points,2), ...
								 't',ctr_t, ...
								 'u',zeros(ctr_points, dcp.m));
% Add the initial condition
ctr_sim.y(1,:) = [(0.1-0.05)*rand(1, 1)+0.05, (0.2-0.05)*rand(1, 1)+0.05];
% Preallocate the vector of inputs
for step = 2 : ctr_points
	% 1. Time vector for the iteration
	t_step = (step-2)*Dt:Dt:(step-2+Np-1)*Dt;
	% 2. Calculate the reference for the Prediction Horizon
	y_ref = [bio_ref(t_step);zeros(1,length(t_step))]'; %
	% 3. Calculate the input for this simulation
	ctr_sim.u(step-1,:) = pqMPCu(ctr_sim.y(step-1,:), ... Output from last iteration
		y_ref, ... reference output for this sim step
		dcp, ... 
		mpc_mat, ...
		[0.01 2] ... vector of delta[u_min u_max]
		);
	% 4. Simulate the system with that input
	[~, y_step] = ode23s(@(t,x)bioODE(t,x, ...
		@(t)ctr_sim.u(step-1)), ... the input just calculated
		t_step(1:2), ... in the first interval of t_step
		ctr_sim.y(step-1,:), ... From the value at the current time
		odeSettings);
	ctr_sim.y(step,:) = y_step(end,:)+ normrnd(0,sigma,size(y_step(end,:)));
	% 5. What happens with the ouput from the EDMD?
end
% Plot the result
%
ctr_fig = figure(3);
clf
tiledlayout(3,1,TileSpacing="tight")
nexttile([2 1])
hold on
outp = plot(ctr_sim.t,ctr_sim.y);
refp = plot(ctr_sim.t,bio_ref(ctr_sim.t));
ylabel('$x$ [g/L]',Interpreter='latex')
legend('Biomass','Substrate','Bio Reference',Interpreter='latex')
title('pqMPC Control',Interpreter='latex')
nexttile
up = plot(ctr_sim.t,ctr_sim.u,Color=[0 .4 .7 trans]);
ylabel('$D$',Interpreter='latex')
xlabel('$t$ [h]', Interpreter='latex')
legend('Dillution rate',Interpreter='latex')
set(gcf,'PaperPosition', [0, 0, 16, 8])
saveas(ctr_fig,"~/Documents/BioReactorConf/figures/ctr_slow_low.fig")
saveas(ctr_fig,"~/Documents/BioReactorConf/figures/ctr_slow_low.eps", "epsc")




