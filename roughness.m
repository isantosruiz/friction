function epsilon = roughness(Re,f)
%ROUGHNESS  Compute the relative roughness coefficient
%   of a pipe from values of the friction factor and
%   Reynolds number for different flow conditions.
%---------------------------------------------------------
%   Sintax
%      epsilon = roughness(Re,f)
%---------------------------------------------------------
%   Arguments
%           Re : Array of Reynolds numbers.
%            f : Array of Darcy-Weisbach friction factors.
%      epsilon : Relative roughness coefficient, k/D.
%---------------------------------------------------------
%   Examples
%      Re = [47525, 74725, 99490, 123013];
%      f = [0.022786, 0.021086, 0.020241, 0.019698]; 
%      epsilon = roughness(Re,f)
%---------------------------------------------------------
%   Ildeberto de los Santos Ruiz, 2018
%   Certified MATLAB Associate
%---------------------------------------------------------
objfunc = @(epsilon)(LambertW(Re(:),epsilon)-f(:))'*(LambertW(Re(:),epsilon)-f(:));
options = optimoptions('fsolve','Algorithm','Levenberg-Marquardt','Display','off');
epsilon = fsolve(objfunc,eps,options);
    function f = LambertW(Re,epsilon)
        a = 2.51./Re;
        b = epsilon/3.7;
        f = 1./(2*lambertw(0,log(10)./(2*a).*10.^(b./(2*a)))/log(10)-b./a).^2;
    end
end