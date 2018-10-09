function [ Wp ] = simpleMeshing( turbData,meshSetup )
    Wp.Lx    = max(turbData.Crx) + meshSetup.distance_N;
    Wp.Ly    = max(turbData.Cry) + meshSetup.distance_E;
    Wp.ldx   = linspace(0,Wp.Lx,meshSetup.Nx);
    Wp.ldy   = linspace(0,Wp.Ly,meshSetup.Ny);
    Wp.ldxx  = repmat(Wp.ldx',1,meshSetup.Ny);
    Wp.ldyy  = repmat(Wp.ldy,meshSetup.Nx,1);
    Wp.ldx2  = 0.5*(Wp.ldx(1:end-1)+Wp.ldx(2:end));
    Wp.ldx2  = [Wp.ldx2 2*Wp.ldx2(end)-Wp.ldx2(end-1)];
    Wp.ldy2  = 0.5*(Wp.ldy(1:end-1)+Wp.ldy(2:end));
    Wp.ldy2  = [Wp.ldy2 2*Wp.ldy2(end)-Wp.ldy2(end-1)];
    Wp.ldxx2 = repmat(Wp.ldx2',1,meshSetup.Ny);
    Wp.ldyy2 = repmat(Wp.ldy2,meshSetup.Nx,1);
    Wp.Nx    = meshSetup.Nx;
    Wp.Ny    = meshSetup.Ny;
end