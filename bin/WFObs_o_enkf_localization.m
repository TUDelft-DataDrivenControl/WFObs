function [ strucObs ] = WFObs_o_enkf_localization( Wp,strucObs,measuredData )
% WFOBS_O_ENKF_LOCALIZATION  Localization function for the EnKF
%
%   SUMMARY
%    This code calculates the physical distance between every state and
%    measurement in the EnKF, and gives a weight to it. Nearby pairs are
%    given a value close to 1, and distant pairs are given a value close to
%    0. This is the localization. Furthermore, this code also scales up
%    these gains through inflation. Localization and inflation are deemed
%    essential for good performance of the EnKF.
%
%   RELEVANT INPUT/OUTPUT VARIABLES
%      see 'WFObs_o.m' for the complete list.
%    

    function kappa = localizationGain( dx, f_locl, l_locl )
        %[ kappa ] = localizationGain( dx, f_locl, l_locl )
        %   This script calculates the localization factor based on the eulerian
        %   distance between two points.
        %
        %    Inputs:
        %     *dx              distance between the two points (m)
        %     *f_locl          localization function ('heaviside' or 'gaspari')
        %     *l_locl          localization length (m)
        %
        %    Outputs:
        %     *kappa           localization factor
        kappa = 0;
        switch lower(f_locl)
            case 'heaviside'
                if dx <= l_locl
                    kappa = 1;
                end;
            case 'gaspari'
                if dx <= l_locl
                    kappa = -(1/4)*(dx/l_locl)^5+(1/2)*(dx/l_locl)^4+...
                        (5/8)*(dx/l_locl)^3 - (5/3)*(dx/l_locl)^2+1;
                elseif dx <= 2*l_locl
                    kappa = (1/12)*(dx/l_locl)^5-(1/2)*(dx/l_locl)^4+(5/8)*(dx/l_locl)^3+...
                        (5/3)*(dx/l_locl)^2-5*(dx/l_locl)+4-(2/3)*(l_locl/dx);
                end;
            otherwise
                error('Wrong localization function specified.');
        end
    end

if strcmp(lower(strucObs.f_locl),'off')
    strucObs.auto_corrfactor = ones(length(measuredData),length(measuredData));
    if strucObs.pe.enabled
        strucObs.cross_corrfactor = ones(strucObs.L,length(measuredData));%/strucObs.M;
    else
        strucObs.cross_corrfactor = ones(strucObs.L,length(measuredData));
    end
else
    disp([datestr(rem(now,1)) ' __  Calculating localization matrices. This may take a while for larger meshes...']);
    rho_locl = struct; % initialize empty structure
    
    % Generate the locations of all model flow states
    stateLocArray = zeros(strucObs.size_output,2);
    for iii = 1:strucObs.size_output
        [~,loci,~]           = WFObs_s_sensors_nr2grid(iii,Wp.mesh);
        stateLocArray(iii,:) = [loci.x, loci.y];
    end
    
    % Generate the locations of all turbines
    turbLocArray = zeros(Wp.turbine.N,2);
    for iii = 1:Wp.turbine.N
        turbLocArray(iii,:) = [Wp.turbine.Crx(iii),Wp.turbine.Cry(iii)];
    end
    
    % Generate the locations of all outputs
    outputLocArray = [];
    for i = 1:length(measuredData)
        if strcmp(measuredData(i).type,'P')
            outputLocArray = [outputLocArray; turbLocArray(measuredData(i).idx,:)];
        elseif strcmp(measuredData(i).type,'u') || strcmp(measuredData(i).type,'v')
            outputLocArray = [outputLocArray; measuredData(i).idx];
        else
            error('You specified an incompatible measurement. Please use types ''u'', ''v'', or ''P'' (capital-sensitive).');
        end
    end
     
%     %%%%% ---- %%%%%%%%%%%%%%%%%%%%    
%     % Covariance localization function for state error covariance matrix
%     rho_locl.autoState = sparse(strucObs.size_output,strucObs.size_output);
%     for iii = 1:strucObs.size_output % Loop over all default states
%         loc1 = stateLocArray(iii,:);
%         for jjj = iii:strucObs.size_output % Loop over all default states
%             loc2 = stateLocArray(jjj,:);
%             dx = sqrt(sum((loc1-loc2).^2)); % displacement between state and output
%             rho_locl.autoState(iii,jjj) = localizationGain( dx, strucObs.f_locl, strucObs.l_locl );
%             rho_locl.autoState(jjj,iii) = rho_locl.autoState(iii,jjj);
%         end
%     end
%     clear iii jjj dx loc1 loc2
%     strucObs.autoState_corrfactor  = rho_locl.autoState;
%     %%%%% ---- %%%%%%%%%%%%%%%%%%%%
    
    % First calculate the cross-correlation between output and states
    if strucObs.se.enabled
        rho_locl.cross = sparse(length(stateLocArray),size(outputLocArray,1));
        for iii = 1:length(stateLocArray) % Loop over all default states
            loc1 = stateLocArray(iii,:);
            for jjj = 1:size(outputLocArray,1) % Loop over all measurements
                loc2 = outputLocArray(jjj,:);
                dx = sqrt(sum((loc1-loc2).^2)); % displacement between state and output
                rho_locl.cross(iii,jjj) = localizationGain( dx, strucObs.f_locl, strucObs.l_locl );
            end
        end
        clear iii jjj dx loc1 loc2
    else
        rho_locl.cross = [];
    end;
    
    % Add cross-correlation between output and model tuning parameters
    if strucObs.pe.enabled
        for iT = 1:length(strucObs.pe.vars)
            if strcmp(strucObs.pe.subStruct{iT},'turbine') % Correlated with all turbines
                crossmat_temp = [];
                for iturb = 1:Wp.turbine.N
                    loc1 = turbLocArray(iturb,:);
                    for jjj = 1:size(outputLocArray,1)
                        loc2 = outputLocArray(jjj,:);
                        dx = sqrt(sum((loc1-loc2).^2)); % displacement between turbine and output
                        crossmat_temp(iturb,jjj) = localizationGain( dx, strucObs.f_locl, strucObs.l_locl );
                    end;
                end;
                if (sum(crossmat_temp,1) <= 0); disp(['Localization too conservative: no correlation between measurements and ' strucObs.pe.vars{iT} '.']); end;
                rho_locl.cross = [rho_locl.cross; max(crossmat_temp)];

            elseif strcmp(strucObs.pe.subStruct{iT},'site') % Correlated with everything in the field equally
%                 rho_locl.cross = [rho_locl.cross; ones(1,strucObs.M)];
                rho_locl.cross = [rho_locl.cross; ones(1,size(outputLocArray,1))/size(outputLocArray,1)]; % Normalized to reduce sensitivity

            else
                disp(['No rules have been set for localization for the online adaption of ' strucObs.pe.vars{iT} '.'])
            end;
        end;
    end
    clear iT dx loc1 loc2 iii jjj crossmat_temp iturb
    
    % Secondly, calculate the autocorrelation of output
    rho_locl.auto = sparse(size(outputLocArray,1),size(outputLocArray,1));
    for iii = 1:size(outputLocArray,1) % for each output
        loc1 = outputLocArray(iii,:);
        for jjj = 1:size(outputLocArray,1) % for each output
            loc2 = outputLocArray(jjj,:);
            dx = sqrt(sum((loc1-loc2).^2));
            rho_locl.auto(iii,jjj) = localizationGain( dx, strucObs.f_locl, strucObs.l_locl );
        end;
    end;
    clear rho_locl_t dx loc1 loc2 iii jjj
    
    % Implement the effect of covariance inflation on localization
    strucObs.auto_corrfactor  = rho_locl.auto;
    strucObs.cross_corrfactor = rho_locl.cross * sqrt(strucObs.r_infl);
end
end