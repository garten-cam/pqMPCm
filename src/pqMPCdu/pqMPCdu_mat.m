classdef pqMPCdu_mat
	%PQMPCMAT recieves a decomposition object, an a set of paramenters to
	%return the necessary matrices for the MPC with integral action that
	%works on an arbitrary pqDecomposition. There are the matrices that are
	%constat throughout the optimization. The remaining matrices depend on
	%the iteration
	
	properties
		Np    % Prediction horizon
		Nc    % Control Horizon
		Atd		% A tilde matrix
		Btd		% B tilde matrix
		Ctd		% C tilde matrix
		Ci		% pseudo inverse of C tilde
		hA		% hat A matrix
		hB		% hat B matrix
		hQx		% hat Qx matrix
		hH		% hat H matrix for the inequalities
		qpH		% H for the quadratic program
		% Aiq % inequality matrix of pq problem
	end
	
	methods
		function obj = pqMPCdu_mat(dcp, Np, Nc, Qy, Qu, Sy)
			%pqMPCmat Calculates the necessary matrices for controlling a
			%dynamical system that has a pqEDMD form.
			% 
			% dcp = decompostion object, from a pqEDMD decomposition
			% Np = prediction Horizon
			% Qx  = 
			obj.Np  = Np; % Save the prediction horizon
			obj.Nc  = Nc; % Save the control horizon
			obj.Atd = [dcp.A, dcp.B; zeros(size(dcp.B')), eye(dcp.m)];
			obj.Btd = [dcp.B; eye(dcp.m)];
			obj.Ctd = [dcp.C, zeros(height(dcp.C),dcp.m)];
			obj.Ci  = pinv(dcp.C);
			% hats
			[obj.hA, obj.hB] = obj.ab_hats(Np, Nc);
			% H matrices
			% obj.hH  = kron(eye(Np), obj.Ctd);
			[obj.hQx, obj.qpH] = obj.quadProgH(dcp, Np, Qy, Qu, Sy);
			obj.hH  = kron(eye(Np), obj.qpH);
			% Inequality 
			% Including x_min and x_max
			% obj.Aineq = [obj.hH*obj.hB;... y_max
			% 	          -obj.hH*obj.hB;... y_min
			% 						 obj.hB       ;... this is for x_max
			% 						-obj.hB];        % this is for x_min
			% Excluding x_min and x_max
			% obj.Aiq = [-obj.hH*obj.hB;... y_max
			% 	          +obj.hH*obj.hB]; % y_min
		end
		function [hQx, qph] = quadProgH(obj, dcp, Np, Qy, Qu, Sy)
			% hat matrices
			Qx = obj.Qy2Qx(dcp,Qy);
			Sx = obj.Qy2Qx(dcp,Sy);
			hQx = blkdiag( ...
				kron(eye(Np-1), blkdiag(dcp.C'*Qx*dcp.C, zeros(dcp.m))), ...
				blkdiag(dcp.C'*Sx*dcp.C, zeros(dcp.m)));
			hQu = kron(eye(Np), Qu);
			qph = obj.hB'*hQx*obj.hB + hQu;
			if ~issymmetric(qph)
				qph = (qph+qph')/2;
			end
		end
		function [ha, hb] = ab_hats(obj, Np, Nc)
			apow = arrayfun(@(Npi){obj.Atd^Npi},0:Np)';
			ha = cell2mat(apow(2:end));
			% Preallocate hat(B)
			hb = arrayfun(@(idx){zeros(size(obj.Btd))},ones(Np));
			% Get the full columns
			bfh = cellfun(@(apwi){apwi*obj.Btd},apow(1:end-1));
			for idx = 1 : Nc
				hb(idx:end,idx) = bfh(1:end-(idx-1));
			end
			hb = cell2mat(hb);
		end
	end
	methods (Static)
		function Xx = Qy2Qx(dcp, Q)
			% obsf = dcp.obs.obs_function;
			% if length(dcp.obs.polynomials_order) < dcp.num_obs
			% 	% We have a bias
			% 	Xx = diag([1 obsf(diag(Q)')]);
			% else
			% 	Xx = diag(obsf(diag(Q)'));
			% end
			Xx = dcp.Cob'*Q*dcp.Cob;
		end
	end
end