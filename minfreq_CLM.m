function [n2, m2, xi2, rhomat, ri, ro, cu, wmin1, wmin2]  = minfreq_CLM(a, z0, H, f, nb,vm, rm, Lr, Lz)

    h = H/2;
    
    v = @(r,z) vm*(r/rm)*exp(-0.5*(r/rm)^a)*exp(-(z-z0)^2/h^2);
    v2 = @(r,z) vm*exp(-0.5*(r/rm)^a)*exp(-(z-z0)^2/h^2);
    zeta = @(r,z) (2-0.5*a*(r/rm)^a)*v(r,z)/r;
    dvz = @(r,z) -2*(z-z0)*v(r,z)/h^2;
    
    Xi2 = @(r,z) (f+2*v(r,z)/r)*(f+2*v(r,z)/r - a*v(r,z)*(r/rm)^a/(2*r));
    M2 = @(r,z) -2*((z-z0)/h^2)*(f+2*v(r,z)/r)*v(r,z);
    
    N2 = @(r,z) 2*v(r,z)*(rm/r)^2/(h^2)*((1-2*(z-z0)^2/h^2)*f*r + (1-4*(z-z0)^2/h^2)*v(r,z)) + nb^2;
    rho = @(r,z) -1000*(2*v2(r,z)*(f*rm + v2(r,z))*z/h^2 + (nb^2)*z)/9.81;
    
    wminp = @(r,z) sqrt((N2(r,z)+Xi2(r,z) + sqrt((N2(r,z)-Xi2(r,z))^2 + 4*M2(r,z)^2))/2);
    wminm = @(r,z) sqrt((N2(r,z)+Xi2(r,z) - sqrt((N2(r,z)-Xi2(r,z))^2 + 4*M2(r,z)^2))/2);
    
    rmat = (1/100:1/100:Lr/rm)*rm; zmat = (-1/100:-1/100:Lz/H)*H;
    
    for i=1:1:length(rmat)
        for j=1:1:length(zmat)
            wmin1(i,j)=wminp(rmat(i),zmat(j));
            wmin2(i,j)=wminm(rmat(i),zmat(j));
            n2(i,j) = N2(rmat(i), zmat(j));
            m2(i,j) = M2(rmat(i), zmat(j));
            xi2(i,j) = Xi2(rmat(i), zmat(j));
            ro(i,j) = zeta(rmat(i), zmat(j))/f;
            ri(i,j) = N2(rmat(i), zmat(j))/abs(dvz(rmat(i), zmat(j)))^2;
            cu(i,j) = 2*v(rmat(i), zmat(j))/(f*rmat(i));
            rhomat(i,j) = rho(rmat(i), zmat(j));
        end
    end
end