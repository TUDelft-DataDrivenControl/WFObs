%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 'WFObs_s_initialize.m'
%  This script loads the model and observer settings. It also prepares
%  the meshing and all variables required for simulation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run(['Configurations/' configName]);% Load configuration file
if strucScript.printProgress
    disp(' WindFarmObserver (WFObs)');
    disp(' ');
    disp(['     Meshing:   ' Wp.name ]);
    disp(['     Observer: ' strucObs.filtertype ]);
    disp(' ');
end;

% Create destination folder & save filter settings
if (strucScript.saveplots + strucScript.saveest + strucScript.saveworkspace > 0)
    mkdir(strucScript.savepath); 
    save([strucScript.savepath '/' strucObs.filtertype '_filtersettings.mat']); 
end; 

% Default settings: following WFSim options are never used in WFObs
options.Projection      = 0;    % Use projection
options.exportLinearSol = 0;    % Export linear solution
options.Derivatives     = 0;    % Calculate derivatives/gradients for system

if strucScript.printProgress
    disp([datestr(rem(now,1)) ' __  Initializing simulation model.']); 
    disp([datestr(rem(now,1)) ' __  Grid-turbine mismatch effects:']); 
    disp(' ')
end;

[Wp,sol,sys,Power,CT,a,Ueffect,input,B1,B2,bc] = ...
    InitWFSim(Wp,options,strucScript.plotMesh); % Initialize model

% Add noise to initial conditions
[sol.u,sol.uu]  = deal(sol.u + randn(Wp.mesh.Nx,Wp.mesh.Ny)*strucObs.noise_init);
[sol.v,sol.vv]  = deal(sol.v + randn(Wp.mesh.Nx,Wp.mesh.Ny)*strucObs.noise_init);

% load a default random seed for consistency
if strucObs.loadrandomseed; load('randomseed'); rng(randomseed); clear randomseed; end;
if strucScript.saveest; save([strucScript.savepath '/' strucObs.filtertype '_est' num2str(datanroffset) '.mat'],'sol'); end; 

% Define what the system should predict (with or without pressures)
strucObs.size_state = Wp.Nu + Wp.Nv + Wp.Np;
if options.exportPressures == 0
    strucObs.size_output = Wp.Nu + Wp.Nv;
else
    strucObs.size_output = Wp.Nu + Wp.Nv + Wp.Np;
end;

% Define measurement locations
sensorsfile = load(sensors_path);
strucObs.obs_array   = unique([sensorsfile.sensors{1}.obsid; sensorsfile.sensors{2}.obsid]);
clear sensorsfile

% Setup blank figure windows
if strucScript.Animation > 0
    if strucScript.plotcontour
        scrsz = get(0,'ScreenSize'); 
        h2=figure('color',[1 1 1],'Position',[50 50 floor(scrsz(3)/1.1) floor(scrsz(4)/1.1)], 'MenuBar','none','ToolBar','none','visible', 'on');
    end; 
    if strucScript.plotpower
        h3=figure;
    end;
    if strucScript.ploterror
        h4=figure;
    end;
end;

% Create global RCM vector
[sysRCM,~,~,~,~] = Make_Ax_b(Wp,sys,sol,input{1},B1,B2,bc,1,options);
sys.pRCM         = sysRCM.pRCM;  clear sysRCM;

klen = length(num2str(Wp.sim.NN));        % used for proper spacing in cmd output window
tlen = length(num2str(Wp.sim.time(end))); % length

if strucScript.printProgress
    disp([datestr(rem(now,1)) ' __  Finished initialization sequence.']);
end;