function [ kappa ] = WFObs_o_enkf_localization( dx, f_locl, l_locl )
%[ kappa ] = WFObs_o_enkf_localization( dx, f_locl, l_locl )
%   This script calculates the localization factor based on the eulerian
%   distance between two points.
%
%    Inputs:
%     *dx              distance between the two points (m)
%     *f_locl          localization function ('heaviside' or 'gaspari')
%     *l_locl          localization length (m)
%
%    Outputs:
%     *kappa           localization factor

kappa = 0;
switch lower(f_locl)
    case 'heaviside'
        if dx <= l_locl
            kappa = 1;
        end;
    case 'gaspari'
        if dx <= l_locl
            kappa = -(1/4)*(dx/l_locl)^5+(1/2)*(dx/l_locl)^4+...
                (5/8)*(dx/l_locl)^3 - (5/3)*(dx/l_locl)^2+1;
        elseif dx <= 2*l_locl
            kappa = (1/12)*(dx/l_locl)^5-(1/2)*(dx/l_locl)^4+(5/8)*(dx/l_locl)^3+...
                (5/3)*(dx/l_locl)^2-5*(dx/l_locl)+4-(2/3)*(l_locl/dx);
        end;
    otherwise
        error('Wrong localization function specified.');
end;
end

