function [sol] = cleanupSolStruct(model,LESData)
    sol = model.sol;
    Wp = model.Wp;
    
    % Use interpolation to determine error between solutions
    LES_sol.u = LESData.uInterpolant(sol.time*ones(size(Wp.mesh.ldxx2(:))),Wp.mesh.ldyy(:),Wp.mesh.ldyy(:));
    LES_sol.v = LESData.vInterpolant(sol.time*ones(size(Wp.mesh.ldxx(:))),Wp.mesh.ldyy2(:),Wp.mesh.ldyy(:));
    
    flowError = [sol.v(:)-LES_sol.v; sol.u(:)-LES_sol.u];
    sol.score = struct('RMSE_flow',rms(flowError),'maxError_flow',max(abs(flowError))); % Errors
    sol.site = model.Wp.site; % Save site info too (contains model parameters that may be estimated)
end

