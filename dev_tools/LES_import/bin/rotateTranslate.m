function [ flowData,turbData,u_Inf,v_Inf,cornerPoints ] = rotateTranslate( flowData,turbData,meshSetup )
% Keep track of corner points
minIdxY = min(flowData.yu)==flowData.yu;
minIdxX = min(flowData.xu)==flowData.xu;
maxIdxY = max(flowData.yu)==flowData.yu;
maxIdxX = max(flowData.xu)==flowData.xu;
blc = [find(minIdxY+minIdxX==2)]; % bottom-left idx
brc = [find(maxIdxY+minIdxX==2)]; % bottom-right idx
trc = [find(maxIdxY+maxIdxX==2)]; % top-right idx
tlc = [find(minIdxY+maxIdxX==2)]; % top-left idx

% Determine freestream conditions flow field by checking all corners
disp('Rotating and translating grid according to u_Inf and v_Inf...')
u_Inf = median(flowData.u(:));
v_Inf = median(flowData.v(:));
WD = atan(v_Inf/u_Inf); % Wind direction [rad]
if abs(WD) > deg2rad(2.5) % Only rotate if mismatch > 2.5 degrees 
    [flowData,turbData] = rotateMesh(flowData,turbData,WD);
    disp('TURBINE YAW ROTATION ONLY VALID FOR (OLD) SOWFA SIMULATIONS!!')
    turbData.phi = -(turbData.phi - 90. + rad2deg(WD)); % rotate to new field (SOWFA only)
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

% Determine corner points of rotated and translated mesh
cornerPoints = [flowData.yu(blc), flowData.xu(blc);...
                flowData.yu(brc), flowData.xu(brc);...
                flowData.yu(trc), flowData.xu(trc);...
                flowData.yu(tlc), flowData.xu(tlc)];
end