function [ grid,loc,element ] = WFObs_s_sensors_nr2grid( stateid, Wp)
%WFObs_s_sensors_nr2grid(stateid,Wp) Converts observation id to grid coordinates.

if stateid <= 0
    error('State id too small. Check consistency with Wp.');
elseif stateid <= (Wp.Nx-3)*(Wp.Ny-2)
    element = 'u';
    grid.x = 3+floor( (stateid-1)/(Wp.Ny-2) );
    grid.y = 1+stateid-(grid.x-3)*(Wp.Ny-2);
    loc.x  = Wp.ldxx2(grid.x,1);
    loc.y  = Wp.ldyy(1,grid.y);
    
elseif stateid <=(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)
    element = 'v';
    stateid  = stateid - (Wp.Nx-3)*(Wp.Ny-2);
    grid.x = 2+floor( (stateid-1)/(Wp.Ny-3));
    grid.y = 2+stateid-(grid.x-2)*(Wp.Ny-3);
    loc.x  = Wp.ldxx(grid.x,1);
    loc.y  = Wp.ldyy2(1,grid.y);
    
elseif stateid <=(Wp.Nx-3)*(Wp.Ny-2)+(Wp.Nx-2)*(Wp.Ny-3)+(Wp.Nx-2)*(Wp.Ny-2)-2
    element = 'p';
    stateid  = stateid - ((Wp.Nx-3)*(Wp.Ny-2) + (Wp.Nx-2)*(Wp.Ny-3));
    grid.x = 2+floor( (stateid-1)/(Wp.Ny-2));
    grid.y = 1+stateid-(grid.x-2)*(Wp.Ny-2);
    loc.x  = Wp.ldxx(grid.x,1);
    loc.y  = Wp.ldyy(1,grid.y);
    
else
    error('State id too large. Check consistency with Wp.');
end;

end

