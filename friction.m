function [f] = friction(q)
    re_friction = q(1)*q(2)*q(3)/q(4);
    Re_c = 3e3;
    
    if (re_friction < Re_c)
        disp(['Laminar, Re = ',num2str(re_friction)]);
        f = 64 / re_friction;
    else
        disp(['Turbulent, Re = ',num2str(re_friction)]);
        eqn = @(x) 1/x +2*log10( (q(5)/3.7/q(3)) + (2.51/re_friction/x));
        f_temp = fsolve(eqn,1e-1,optimoptions('fsolve','Display','off'));
%         disp(['sqrt(f) = ',num2str(f_temp)]);
        f = f_temp^2;
    end

end