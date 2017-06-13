function [ measured ] = WFObs_s_preprocess( Wp, measured )
%% Preprocess measurements to determine freestream wind direction and speed

    function [alpha] = getAlpha(x0,y0,x1,y1,D)
        dx = abs(x1-x0);
        dy = abs(y1-y0);
        L = sqrt(dx*dx+dy*dy);
        alpha = 1.3*atan(2.5*D/L+0.15) + pi/180. * 10.;
        alpha = alpha*180./pi
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
            gamma = abs(beta - 360. - theta) - alpha/2.
        else
            gamma = np.abs(beta - theta) - alpha/2.
        end;
    end
        
    function [isWaked] = isWaked(x0,y1,x1,y1,D,wd)
        alpha = getAlpha(x0,y0,x1,y1,D);
        beta  = getBeta(x0,y0,x1,y1);
        gamma = getGamma(beta,wd,alpha);
        if gamma < 0
            isWaked = 1;
        else
            isWaked = 0;
        end;
    end

    function [isWakedAny] = isWakedAny(turbineIdx,x0,y1,x1,y1,D,wd)
        alpha = getAlpha(x0,y0,x1,y1,D);
        beta  = getBeta(x0,y0,x1,y1);
        gamma = getGamma(beta,wd,alpha);
        if gamma < 0
            isWaked = 1;
        else
            isWaked = 0;
        end;
    end

end