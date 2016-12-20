function [ sol,strucObs ] = WFObs_o_exkf(strucObs,Wp,sys,B1,B2,bc,input,measured,sol,k,it,options)
%[ sol, strucObs ] = WFObs_o_exkf(strucObs,Wp,sys,B1,B2,bc,input,measured,sol,k,it,startUniform)
%   This script calculates the optimally estimated system state vector
%   according to the measurement data and the internal model using the
%   Extended Kalman Filter (ExKF).
%
%    Inputs:
%     *strucObs        structure containing observer information for time k-1
%     *Wp              structure containing meshing information
%     *sys             structure containing system matrices
%     *B1,B2,bc        system matrices related to the boundary conditions
%     *input           structure with turbine control settings (yaw and a)
%     *measured        structure containing (SOWFA) measurement information
%     *sol             flow fields and system state vector for time k-1
%     *k,it            current sample and iteration number, respectively
%     *options         structure containing model/script option information
%
%    Outputs:
%     *sol             flow fields and system state vector for time k
%     *strucObs        structure containing observer information for time k

    if k == 1
        % Setup covariance and system output matrices
        if options.exportPressures
            strucObs.Pk    = sparse(eye(strucObs.size_state))*strucObs.P_0;
            strucObs.Htt   = sparse(eye(strucObs.size_state));
            strucObs.Htt   = strucObs.Htt(strucObs.obs_array,:);
            strucObs.Q_k   = strucObs.Q_k*eye(strucObs.size_state);
        else
            strucObs.Pk    = sparse(eye(strucObs.size_output))*strucObs.P_0;
            strucObs.Htt   = sparse(eye(strucObs.size_output));
            strucObs.Htt   = strucObs.Htt(strucObs.obs_array,:);
            strucObs.Q_k   = strucObs.Q_k*eye(strucObs.size_output);
        end;
        
        % Calculate the sparsification matrix for system matrix F
        if strucObs.sparseF
            [sysf,power,~,~,~]    = Make_Ax_b(Wp,sys,sol,input,B1,B2,bc,k,options); % Create system matrices
            Fk(sys.pRCM,sys.pRCM) = sysf.A(sys.pRCM,sys.pRCM)\sysf.Al(sys.pRCM,sys.pRCM); % Linearized A-matrix at time k
            strucObs.indFsparse   = abs(Fk)>1e-2;
            if ~options.exportPressures
                strucObs.indFsparse = strucObs.indFsparse(1:strucObs.size_output,1:strucObs.size_output);
            end;
        end;
    end;
    
    % Calculate forecasted state vector
    [sysf,power,~,~,~]    = Make_Ax_b(Wp,sys,sol,input,B1,B2,bc,k,options);       % Create system matrices
    Fk(sys.pRCM,sys.pRCM) = sysf.A(sys.pRCM,sys.pRCM)\sysf.Al(sys.pRCM,sys.pRCM); % Linearized A-matrix at time k
    Bk(sys.pRCM,:)        = sysf.A(sys.pRCM,sys.pRCM)\sysf.Bl(sys.pRCM,:);        % Linearized B-matrix at time k    
        
    if ~options.exportPressures % Neglect pressure terms
        Fk = Fk(1:strucObs.size_output,1:strucObs.size_output);
        Bk = Bk(1:strucObs.size_output,:);                     
    end;
    
    if strucObs.sparseF % Enforce sparsification
        Fk = strucObs.indFsparse .* Fk; 
    end;
    
    [solf,~] = Computesol(sysf,input,sol,k,it,options);  % Compute forecasted solution by standard WFSim function
    Pf       = Fk*strucObs.Pk*Fk' + strucObs.Q_k;        % Covariance matrix P for x(k) knowing y(k-1)
    if strucObs.diagP; Pf = diag(diag(Pf)); end;         % Enforce sparsification of Pf
    
    if ~options.exportPressures
        solf.x = solf.x(1:strucObs.size_output);              % Remove pressure terms
    end;
    
    % Analysis update
    Kgain       = Pf(:,strucObs.obs_array)*pinv(Pf(strucObs.obs_array,strucObs.obs_array)+strucObs.R_k); % Kalman gain
    sol.x       = solf.x + Kgain*(measured.sol(strucObs.obs_array)-solf.x(strucObs.obs_array));          % Optimally predicted state vector
    strucObs.Pk = (eye(size(Pf))-Kgain*strucObs.Htt)*Pf;                                                 % State covariance matrix
    if strucObs.diagP; strucObs.Pk = diag(diag(strucObs.Pk)); end;  % Enforce sparsification of Pk
end

