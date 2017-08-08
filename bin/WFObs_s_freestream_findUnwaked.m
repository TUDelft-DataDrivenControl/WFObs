function [ unwakedTurbines ] = WFObs_s_determineUpstream( Wp, wd )
% WFOBS_S_DETERMINEUPSTREAM  Determine which turbines are upstream
%
%   SUMMARY
%    This code uses a very simple static wake model to determine
%    approximately the width of wakes. This information, in combination
%    with the freestream wind direction, allows one to determine which
%    turbines are operating in the freestream. Credits go to Paul Fleming
%    from NREL, who took it from another paper (..?)
%
%   RELEVANT INPUT/OUTPUT VARIABLES
%     - Wp: this struct contains all the simulation settings related to the
%           wind farm, the turbine inputs, the atmospheric properties, etc.
%
%     - Wd: the freestream wind direction, estimated using wind vane
%           measurements.
%

% Define necessary subfunctions
    function [alpha] = getAlpha(x0,y0,x1,y1,D)
        dx = abs(x1-x0);
        dy = abs(y1-y0);
        L = sqrt(dx*dx+dy*dy);
        alpha = 1.3*atan(2.5*D/L+0.15) + pi/180. * 10.;
        alpha = alpha*180./pi;
    end

    function [beta] = getBeta(x0,y0,x1,y1)
        beta = atan2(x0-x1,y0-y1)*180./pi;
        if beta < 0.
            beta = beta + 360.;
        end;
    end

    function [gamma] = getGamma(beta,theta,alpha)
        if((beta >= 0. && beta < 90.) && (theta > 270. && theta < 360.))
            gamma = abs(beta + 360. - theta) - alpha/2.;
        elseif ((beta > 270. && beta < 360.) && (theta >= 0. && theta < 90.))
            gamma = abs(beta - 360. - theta) - alpha/2.;
        else
            gamma = abs(beta - theta) - alpha/2.;
        end;
    end

    function [isWaked] = isWaked(x0,y0,x1,y1,D,wd)
        alpha = getAlpha(x0,y0,x1,y1,D);
        beta  = getBeta(x0,y0,x1,y1);
        gamma = getGamma(beta,wd,alpha);
        if gamma < 0
            isWaked = 1;
        else
            isWaked = 0;
        end;
    end

    function [isWakedAny] = isWakedAny(turbineIdx,D,wd,x,y)
        x1 = x(turbineIdx);
        y1 = y(turbineIdx);
        waked = 0;
        
        % Check if waked by any turbine
        for i = 1:length(x)
            if i ~= turbineIdx
                waked = waked + isWaked(x(i),y(i),x1,y1,D,wd);
            end;
        end;
        
        % If waked by any turbine, return True
        isWakedAny = (waked > 0.5);
    end

% Determine unwaked turbines
unwakedTurbines = [];
x  = Wp.turbine.Crx;
y  = Wp.turbine.Cry;
D  = Wp.turbine.Drotor;

for j = 1:length(x)
    if isWakedAny(j,D,wd,x,y) == false
        unwakedTurbines = [unwakedTurbines, j];
    end
end
end