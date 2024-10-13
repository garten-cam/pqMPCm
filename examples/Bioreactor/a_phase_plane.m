% Script for the simulation of the simple bioreactor
% simulate the phase plane and the approximation by the pqEDMD
rng(1234)
% define the parameters for the simulation
num_ics = 60; % Number of initial conditions for the test
% Create the initial conditions for the samples
% Even though the initial conditions should always start at the lower left
% corner, this will plot the phase plane with constant dilution rate.
% (b-a)*rand(n,m)+a
% For $x_0\in (0,\;1.58]\times(0,\;4]$
ics = [
	(1.58-0.01)*rand(num_ics, 1)+0.01, (4-0.01)*rand(num_ics, 1)+0.01
	];
%
tfin = 40;
n_points = 1*tfin+1;
%
% The first figure is the phase plane with dilution rate D=O.2
in = 0.2*ones(num_ics,1);
% Preallocate the experiments if the samples are not in a column
% vector, the algorithm breaks. Dammit
exp = arrayfun(@(z)struct('y', zeros(n_points, 2), ...
	'u', zeros(n_points, 1),... forcing signals
	't', zeros(n_points, 1)),1:num_ics)';


odeSettings = odeset('RelTol',1e-3,'AbsTol',1e-6);
for orb = 1 : num_ics
	[exp(orb).t, exp(orb).y] = ode23s(@(t,x)cstr_ode(t,x, ...
		in(orb)),...
		linspace(0,tfin,n_points), ...
		ics(orb,:), ...
		odeSettings);
	exp(orb).u = in(orb)*ones(size(exp(orb).t));
end

%
ppfig = figure(1);
clf
hold on
ppp = arrayfun(@(ex)plot(ex.y(:,1),ex.y(:,2),'b'),exp);
xlabel('$x_1$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
title("Phase plane", Interpreter="latex")
set(gcf,'PaperPosition', [0, 0, 7, 7])

saveas(ppfig,"~/Documents/BioReactorConf/figures/a_phase_plane.fig")
saveas(ppfig,"~/Documents/BioReactorConf/figures/a_phase_plane.eps", "epsc")
