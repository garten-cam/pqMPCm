% function u = pqMPCu(y, y_ref, dcp, mpc_mat, y_mx, du_mx)
function u = pqMPCu(y, y_ref, dcp, mpc_mat, u_mx)
%PQMPCU Calculates the optimin input u for the step that brings y to its
%desired setpoint
%
% arguments:
% y				: output from the last step
% y_ref		: reference that the output must follow
% dcp			: decomposition object
% mpc_mat : MPC matrices object
% u_mx		: limits for the input, min and max bounds.
% Calculate the current state
x = mpc_mat.Ci*vlift(dcp, y)';
% Calculate the setpoints of the states
xr = xref(dcp, mpc_mat, y_ref);
% Deviation of the state
x_ev = mpc_mat.hA*x; % Evolution of the current state for Np steps
xd = x_ev - xr;
% h part of the quadratic program u'Hu + h'u
h = xd'*mpc_mat.hQx*mpc_mat.hB;
% Inequalities
% biq = [repmat(vlift(dcp, y_mx(:,2)')',mpc_mat.Np,1)-mpc_mat.hH*x_ev;...
	    % -repmat(vlift(dcp, y_mx(:,1)')',mpc_mat.Np,1)+mpc_mat.hH*x_ev;...
		  % ];
% Upper and lower bounds
ub = repmat(u_mx(:,2),mpc_mat.Np,1);
lb = repmat(u_mx(:,1),mpc_mat.Np,1);
% Ready... All the igredients are ready
% Options
alg_options = optimoptions('quadprog','Algorithm','interior-point-convex', Display='off'); 
[uv, ~, ex_fl] = quadprog(mpc_mat.qpH, ... H
	                        h', ...           h matlab calls it f
													[], [], ...     A and B inequality 
													[], [], ...     A and B equality
													lb, ub, ...     Lower and upper bounds
													[], alg_options); 
u = uv(1,:);
% Done, check for flags and finish
if ex_fl < 0
	disp('infeasible qp problem')
end
end
% Some auxiliary functions. this is not a class because it just spits a 'u'
% an 'u'? no idea, I do not need to istantiate an object. I just do not
% want a lot of clutter in a sigle function
function vl = vlift(dcp, v)
% lift the variable v into the function space
obsf = dcp.obs.obs_function;
if length(dcp.obs.polynomials_order) < dcp.num_obs
	vl = [ones(height(v),1) obsf(v)];
else
	vl = obsf(v);
end
end
function xr = xref(dcp, mpc_mat, y_ref)
ylift = vlift(dcp, y_ref);
xr = mpc_mat.Ci*ylift';
xr = reshape(xr,[],1);
end