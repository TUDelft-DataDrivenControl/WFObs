function [ obs_ids ] = WFObs_s_sensors_grid2nr( gridcor, Wp, winddir)
%WFObs_s_sensors_grid2nr(gridcor,Wp,winddir) Converts grid coordinates to obsv vector entry numbers.

switch winddir % Feasibility check and determine observability vector
    case {'u','U'}
        if ( (sum(gridcor(:,1)<=2) > 0) | (sum(gridcor(:,1)>=Wp.mesh.Nx) > 0) ...
                | (sum(gridcor(:,2)<=1) > 0) | (sum(gridcor(:,2)>=Wp.mesh.Ny) > 0) )
            errordlg('Error: you specified sensors on BC grid points. Please retry with different sensor locations.');
            error('Error: you specified sensors on BC grid points. Please retry with different sensor locations.');
        end;
        for j = 1:size(gridcor,1)
            obs_ids(j,1) = (gridcor(j,1)-3)*(Wp.mesh.Ny-2)+gridcor(j,2)-1;
        end;
        
    case {'v','V'}
        if ( (sum(gridcor(:,1)<=1) > 0) | (sum(gridcor(:,1)>=Wp.mesh.Nx) > 0) ...
                | (sum(gridcor(:,2)<=2) > 0) | (sum(gridcor(:,2)>=Wp.mesh.Ny) > 0) )
            errordlg('Error: you specified sensors on BC grid points. Please retry with different sensor locations.');
            error('Error: you specified sensors on BC grid points. Please retry with different sensor locations.');
        end;
        for j = 1:size(gridcor,1)
            obs_ids(j,1) = (gridcor(j,1)-2)*(Wp.mesh.Ny-3)+(gridcor(j,2)-2) + (Wp.mesh.Ny-2)*(Wp.mesh.Nx-3);
        end;
        
    otherwise
        error('Wrong wind direction specified. Please define as "u" or "v".');
end;
end

