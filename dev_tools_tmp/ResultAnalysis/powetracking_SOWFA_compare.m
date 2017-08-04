clear all; clc;

addpath ..\..\export_fig;
addpath ..\SOWFA; % to run ImportSUPERCONOUT

% No turbulence
SCO1 = importSuperCONOUT('M:\3me\dcsc\DataDriven\Data\SOWFA\Yawexcitationcase1\superCONOUT.csv');
SCO2 = importSuperCONOUT('M:\3me\dcsc\DataDriven\Data\SOWFA\Yawexcitationcase2\superCONOUT.csv')
% Turbulence
SCO3 = importSuperCONOUT('M:\3me\dcsc\DataDriven\Data\SOWFA\Yawexcitationcase3\superCONOUT.csv');

P1 = [SCO1.data{1}(:,2) SCO1.data{2}(:,2)];
P2 = [SCO2.data{1}(:,2) SCO2.data{2}(:,2)];
P3 = [SCO3.data{1}(:,2) SCO3.data{2}(:,2)];

filt  = tf([(0.04*2*pi)^2],[1 2*(0.04*2*pi)*0.7 (0.04*2*pi)^2]);
filtd = c2d(filt,1);
P1    = lsim([filtd 0; 0 filtd],SCO1.data{1}(:,2);