function [ flowDataOut,turbDataOut ] = rotateMesh( flowDataIn,turbDataIn,WD )

% Rotation in the frame from x to x', positive
%
%   / \     _ x'
%    |      /|
% x  | WD  /
%    |    /
%    |   /
%    |  /
%    |_/_ _ _ _ > y
%     \
%     _\| y'
%

flowDataOut = flowDataIn;
turbDataOut = turbDataIn;

% Rotated turbine locations
turbDataOut.Crx = cos(WD)*turbDataIn.Crx+sin(WD)*turbDataIn.Cry;
turbDataOut.Cry = cos(WD)*turbDataIn.Cry-sin(WD)*turbDataIn.Crx;

% Rotated mesh
flowDataOut.xu = +cos(WD) * flowDataIn.xu + sin(WD) * flowDataIn.yu;
flowDataOut.yu = -sin(WD) * flowDataIn.xu + cos(WD) * flowDataIn.yu;
flowDataOut.xv = +cos(WD) * flowDataIn.xv + sin(WD) * flowDataIn.yv;
flowDataOut.yv = -sin(WD) * flowDataIn.xv + cos(WD) * flowDataIn.yv;
% figure;plot(flowDataOut.yu,flowDataOut.xu,'.'); hold on; plot(turbDataOut.Cry,turbDataOut.Crx,'ro')
end