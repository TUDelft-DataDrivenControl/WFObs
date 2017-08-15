function [ flowData,turbData,u_Inf,v_Inf ] = rotateTranslate( flowData,turbData,meshSetup )
% Determine freestream conditions flow field by checking all corners
disp('Rotating and translating grid according to u_Inf and v_Inf...')
u_Inf = median(flowData.u(:));
v_Inf = median(flowData.v(:));
WD = atan(v_Inf/u_Inf); % Wind direction [rad]
if abs(WD) > deg2rad(2.5) % Only rotate if mismatch > 2.5 degrees 
    [flowData,turbData] = rotateMesh(flowData,turbData,WD);
    turbData.phi = turbData.phi*pi/180; % from deg to radians
    turbData.phi = turbData.phi - pi/2 + WD; % rotate to new field
end

% Translation
[~,UpstrIndx] = min(turbData.Crx); % Find most upstream turbine
flowData.xv = flowData.xv   - turbData.Crx(UpstrIndx) + meshSetup.distance_S;
flowData.xu = flowData.xu   - turbData.Crx(UpstrIndx) + meshSetup.distance_S;
flowData.yu = flowData.yu   - turbData.Cry(UpstrIndx) + meshSetup.distance_W;
flowData.yv = flowData.yv   - turbData.Cry(UpstrIndx) + meshSetup.distance_W;
turbData.Crx = turbData.Crx - turbData.Crx(UpstrIndx) + meshSetup.distance_S;
turbData.Cry = turbData.Cry - turbData.Cry(UpstrIndx) + meshSetup.distance_W;

% Align freestream conditions 
u_Inf = sqrt(u_Inf^2 + v_Inf^2);
v_Inf = 0;
end