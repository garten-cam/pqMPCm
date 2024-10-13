classdef pqMPCu_mat
	%PQMPCMAT recieves a decomposition object, an a set of paramenters to
	%return the necessary matrices for the MPC with integral action that
	%works on an arbitrary pqDecomposition. There are the matrices that are
	%constat throughout the optimization. The remaining matrices depend on
	%the iteration
	
	properties
		Np    % Prediction horizon
		Nc    % Control Horizon
		Ci		% pseudo inverse of C tilde
		hA		% hat A matrix
		hB		% hat B matrix
		hQx		% hat Qx matrix
		hH		% hat H matrix for the inequalities
		qpH		% H for the quadratic program
		% Aiq % inequality matrix of pq problem
	end
	
	methods
		function obj = pqMPCu_mat(dcp, Np, Nc, Qy, Qu, Sy, x_method)
			%pqMPCmat Calculates the necessary matrices for controlling a
			%dynamical system that has a pqEDMD form.
			% 
			% dcp = decompostion object, from a pqEDMD decomposition
			% Np = prediction Horizon
			% Qc = control Horizon
			% Qy = output
			obj.Np  = Np; % Save the prediction horizon
			obj.Nc  = Nc; % Save the cotrol horizon
			obj.Ci  = pinv(dcp.C);
			% hats
			[obj.hA, obj.hB] = obj.ab_hats(dcp, Np, Nc);
			% H matrices
			obj.hH  = kron(eye(Np), dcp.C);
			[obj.hQx, obj.qpH] = obj.quadProgH(dcp, Np, Qy, Qu, Sy, x_method);
			% Inequality 
			% Including x_min and x_max
			% obj.Aineq = [obj.hH*obj.hB;... y_max
			% 	          -obj.hH*obj.hB;... y_min
			% 						 obj.hB       ;... this is for x_max
			% 						-obj.hB];        % this is for x_min
			% Excluding x_min and x_max
			% obj.Aiq = [obj.hH*obj.hB;... y_max
				          % -obj.hH*obj.hB]; % y_min
		end
		function [hQx, qph] = quadProgH(obj, dcp, Np, Qy, Qu, Sy, x_method)
			% hat matrices
			Qx = obj.Qy2Qx(dcp, Qy, x_method);
			Sx = obj.Qy2Qx(dcp, Sy, x_method);
			hQx = blkdiag(kron(eye(Np-1), dcp.C'*Qx*dcp.C), dcp.C'*Sx*dcp.C);
			hQu = kron(eye(Np), Qu);
			% qph = round((obj.hB'*hQx*obj.hB + hQu)*10^9)/10^9;
			qph = obj.hB'*hQx*obj.hB + hQu;
			if ~issymmetric(qph)
				qph = (qph+qph')/2;
			end
		end
	end
	methods (Static)
		function [ha, hb] = ab_hats(dcp, Np, Nc)
			apow = arrayfun(@(Npi){dcp.A^Npi},0:Np)';
			ha = cell2mat(apow(2:end));
			% Preallocate hat(B)
			hb = arrayfun(@(idx){zeros(size(dcp.B))},ones(Np));
			% Get the full columns
			bfh = cellfun(@(apwi){apwi*dcp.B},apow(1:end-1));
			for idx = 1 : Nc
				hb(idx:end,idx) = bfh(1:end-(idx-1));
			end
			hb = cell2mat(hb);
		end
		function Xx = Qy2Qx(dcp, Q, x_method)
			switch x_method
				case "Cob"
					Xx = dcp.Cob'*Q*dcp.Cob;
				case "Obs"
					obsf = dcp.obs.obs_function;
					if length(dcp.obs.polynomials_order) < dcp.num_obs
						% We have a bias
						Xx = diag([1 obsf(diag(Q)')]);
					else
						Xx = diag(obsf(diag(Q)'));
					end
			end
		end
	end
end