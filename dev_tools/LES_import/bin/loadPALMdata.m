function [ flowDataRaw,turbDataRaw ] = loadPALMdata( filesInFolder,hubHeight )
% Internal function for loading turbine data
    function turbDataOut = loadTurbData(fileName)
        M = importdata(fileName);
        if class(M) == 'struct'
            M = M.data; % Compatibility with different headers
        end
        
        % [Time   UR  Uinf  Ct_adm  a Yaw Thrust Power  WFPower]
        turbDataOut.time  = M(:,1)-M(1,1); % Normalized to start at 0
        turbDataOut.Ur    = M(:,2);
        turbDataOut.Uinf  = M(:,3);
        turbDataOut.CT    = M(:,4);
        turbDataOut.a     = M(:,5);
        turbDataOut.phi   = M(:,6);
        turbDataOut.power = M(:,8);
    end

% Internal function for loading flow data
    function flowDataOut = loadFlowData(fileName)
        flowDataOut.time  = double(nc_varget(fileName,'time'));
        flowDataOut.time  = flowDataOut.time-flowDataOut.time(1); % Normalize to start at 0
        flowDataOut.u     = double(nc_varget(fileName,'u'));
        flowDataOut.v     = double(nc_varget(fileName,'v'));
        flowDataOut.xu    = double(nc_varget(fileName,'xu'));
        flowDataOut.xv    = double(nc_varget(fileName,'x'));
        flowDataOut.yu    = double(nc_varget(fileName,'y'));
        flowDataOut.yv    = double(nc_varget(fileName,'yv'));
        % flowDataOut.zw_3d = double(nc_varget(fileName,'zw_3d'));
        flowDataOut.zu_3d = double(nc_varget(fileName,'zu_3d'));
    end

addpath(genpath('bin/mexcdf'));

% Loop over all files
for jFile = filesInFolder
    jFile = jFile{1};
    
    % Check if it is a flow data
    if ~isempty(strfind(jFile,'.nc')) % if contains .nc
        if exist('flowDataRaw')
            error('Confused by the existence of multiple .nc files.');
        else
            flowDataRaw = loadFlowData(jFile);
            
        end
        % Otherwise, check if it is a turbine file
    elseif ~isempty(strfind(jFile,'turbine_param')) % if contains 'turbine_param'
        stringIdx = strfind(jFile,'turbine_param');
        turbineId = regexp(jFile(stringIdx:end),'\d*','Match');
        turbineId = str2double(turbineId); % Convert to turbine Id
        
        % Import from external file
        turbDataRaw(turbineId) = loadTurbData(jFile);
    else
        disp(['Skipping file with name ' jFile '.'])
    end
end

% Convert array of structs 'turbData' to single struct 'turbData'
fieldNamesListTurb = fieldnames(turbDataRaw);
for j = 1:length(fieldNamesListTurb)
    jName = fieldNamesListTurb{j};
    turbDataTmp.(jName) = [turbDataRaw.(jName)];
end
turbDataRaw      = turbDataTmp;
turbDataRaw.time = turbDataRaw.time(:,1); % Make vector

% Find z-index closest to turbine hub-height to export horizontal slice
nz = round(interp1(flowDataRaw.zu_3d,1:length(flowDataRaw.zu_3d),hubHeight));
flowDataRaw.zu_3d = flowDataRaw.zu_3d(nz);

% Export horizontal slice and make u,v a 3D tensor: (k,y,z)
flowDataRaw.u = squeeze(flowDataRaw.u(:,nz,:,:));
flowDataRaw.v = squeeze(flowDataRaw.v(:,nz,:,:));

% Vectorize everything to support any type of grid
[yu_mesh,xu_mesh] = meshgrid(flowDataRaw.yu,flowDataRaw.xu);
[yv_mesh,xv_mesh] = meshgrid(flowDataRaw.yv,flowDataRaw.xv);
NN     = size(flowDataRaw.u,1);
Nu     = size(flowDataRaw.u,2)*size(flowDataRaw.u,3);
yu_vec = reshape(yu_mesh,Nu,1);
xu_vec = reshape(xu_mesh,Nu,1);
yv_vec = reshape(yv_mesh,Nu,1);
xv_vec = reshape(xv_mesh,Nu,1);

[uvec,vvec] = deal(zeros(NN,Nu));
for j = 1:NN
    uvec(j,:) = reshape(squeeze(flowDataRaw.u(j,:,:))',Nu,1); % vectorized
    vvec(j,:) = reshape(squeeze(flowDataRaw.v(j,:,:))',Nu,1); % vectorized
end
flowDataRaw.u = uvec;
flowDataRaw.v = vvec;
flowDataRaw.xv = xv_vec;
flowDataRaw.xu = xu_vec;
flowDataRaw.yv = yv_vec;
flowDataRaw.yu = yu_vec;
end