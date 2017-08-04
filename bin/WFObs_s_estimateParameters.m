    %         %% Pre-processing: parameter estimation
    %         paramUpdateFreq = 200;
    %         if ~rem(k,paramUpdateFreq)
    %             for i = 1:paramUpdateFreq
    %                 t_tmp               = Wp.sim.time(i+k-paramUpdateFreq+1);
    %                 measured_tmp        = WFObs_s_loadmeasurements(sourcepath,datanroffset,t_tmp,0);
    %                 solq_tmp(:,i)       = measured.solq(strucObs.obs_array);
    %                 solq_tAvg           = mean(solq_tmp,2); % Column average
    %                 
    %                 input_tmp.beta(:,i) = [input{i+k-paramUpdateFreq}.beta];
    %                 input_tmp.phi(:,i)  = [input{i+k-paramUpdateFreq}.phi];
    %                 input_tAvg.beta     = mean( input_tmp.beta,2 );
    %                 input_tAvg.phi      = mean( input_tmp.phi, 2 );
    %                 
    %                 Wp_tAvg             = Wp;
    %                 Wp_tAvg.site.u_Inf  = mean(Wp.saved.u_Inf(k-paramUpdateFreq:k));
    %                 Wp_tAvg.site.v_Inf  = mean(Wp.saved.v_Inf(k-paramUpdateFreq:k));
    %                 Wp_tAvg.sim.h       = Inf; % Steady-state simulation
    %                 
    %                 % Apply changed boundary conditions to update system matrices
    %                 [B1_tAvg,B2_tAvg,bc_tAvg] = Compute_B1_B2_bc(Wp_tAvg);
    %                 B2_tAvg                   = 2*B2_tAvg;
    %             end
    %             clear i sol_tmp t_tmp measured_tmp % Remove old variables
    %             
    %             % Cost function
    %             function J = costFunction(x,Wp_tAvg,input_tAvg,solq_tAvg,B1_tAvg,B2_tAvg,bc_tAvg )
    %                 global strucObs sys options
    %                 Wp_tAvg.turbine.forcescale = x(1);
    %                 
    %                 it        = 0;
    %                 eps       = 1e19;
    %                 epss      = 1e20;
    % 
    %                 while ( eps>conv_eps && it<max_it && eps<epss )
    %                     it   = it+1;
    %                     epss = eps;
    %                     
    %                     %% Calculate optimal solution according to filter of choice
    %                     [sol,Wp,strucObs]   = WFObs_o(strucObs,Wp_tAvg,sys,B1_tAvg,B2_tAvg,bc_tAvg,input_tAvg,[],struct(),1,it_tmp,options);  % Perform the observer update
    %                     [sol,eps]           = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);   % Map solution to flowfields
    %                 end;
    %             end
    %             
    %         end