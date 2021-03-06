function [sol_out] = formatSol(model,sol_true)
    sol_in = model.sol;
    Wp = model.Wp;
    
    flowError = [sol_in.v(:)-sol_true.v(:); sol_in.u(:)-sol_true.u(:)];
    
    sol_out = struct(...
        'k',sol_in.k,...
        'time',sol_in.time,...
        'uEst',sol_in.u,...
        'uTrue',reshape(sol_true.u,Wp.mesh.Nx,Wp.mesh.Ny),...
        'vEst',sol_in.v,...
        'vTrue',reshape(sol_true.v,Wp.mesh.Nx,Wp.mesh.Ny),...
        'pEst',sol_in.p,...
        'PEst',sol_in.turbine.power,...
        'PTrue',sol_true.P,...
        'turbInput',sol_in.turbInput,...
        'measuredData',sol_in.measuredData,...
        'site',Wp.site,...
        'RMSE_flow',rms(flowError),...
        'maxError_flow',max(abs(flowError)),...
        'CPUtime',sol_in.CPUtime...
        );
end

