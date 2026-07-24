%Function for ray tracing but for CLM vortex 
function [r_ray,z_ray,l_ray,m_ray,tmat,omegamat,cgr_ray,cgz_ray,chmat] = raytrace_CLM_ode23(rstart,zstart,t,Lr,Lz,omega0,m0,ch,a,vm,rm,z0,h,f,nb)
    
    %r,z - radial and axial coordinates of the ray
    %l,m - radial and axial wavenumbers
    %cgr,cgz - group velocities in radial and axial velocities
    
    %defining matrices and values
    tmat = []; r_ray = []; z_ray = []; l_ray = []; m_ray = []; chmat = []; 
    alpha_ray = []; cgr_ray = []; cgz_ray = []; omegamat = []; err = 10^(-8);
    ar1mat = []; ar2mat = [];

    % base flow quantities
    v = @(r,z) vm*(r/rm)*exp(-0.5*(r/rm)^a)*exp(-(z-z0)^2/h^2);
    xi2 = @(r,z) (f+2*v(r,z)/r)*(f+2*v(r,z)/r - a*v(r,z)*(r/rm)^a/(2*r));
    M2 = @(r,z) -2*((z-z0)/h^2)*(f+2*v(r,z)/r)*v(r,z);
    n2 = @(r,z) 2*v(r,z)*(rm/r)^2/(h^2)*((1-2*(z-z0)^2/h^2)*f*r + (1-4*(z-z0)^2/h^2)*v(r,z)) + nb^2;

    %gradients of base flow quantities
    dxi2r = @(r,z) (v(r,z)/rm^2)*(f*(r^2/rm^2 - 6) + 4*(v(r,z)/r)*(r^2/rm^2 - 3));
    dxi2z = @(r,z) (2*v(r,z)/r)*((z-z0)/h^2)*(f*(r^2/rm^2 - 4) + 4*(v(r,z)/r)*(r^2/rm^2 - 2));

    dM2r = @(r,z) 2*((z-z0)/h^2)*(v(r,z)/r)*(f*((r/rm)^2 - 1) + 4*v(r,z)*r/rm^2 - 2*v(r,z)/r);
    dM2z = @(r,z) -2*(v(r,z)/r )*(2*(v(r,z)/h^2)*(1-4*(z-z0)^2/h^2) + f*r*(1-2*(z-z0)^2/h^2)/h^2);

    dn2r = @(r,z) -2*(v(r,z)/r)*(2*(v(r,z)/h^2)*(1-4*(z-z0)^2/h^2) + f*r*(1-2*(z-z0)^2/h^2)/h^2);
    dn2z = @(r,z) 4*(rm/r)^2*((z-z0)/h^2)*v(r,z)*(f*r*(-3+2*(z-z0)^2/h^2)/h^2 + 2*v(r,z)*(-3+4*(z-z0)^2/h^2)/h^2);

    %frequency, aspect ratio, group velocity and gradient of frequency
    omega = @(r,z,l,m) (m/sqrt(m^2+l^2))*sqrt(n2(r,z)*(l^2)/(m^2) - 2*M2(r,z)*l/m + xi2(r,z));
    ar1 = @(r,z) (2*M2(r,z) + sqrt(4*M2(r,z)^2 - 4*(n2(r,z)-omega0^2)*(xi2(r,z)-omega0^2)))/(2*(n2(r,z)-omega0^2));
    ar2 = @(r,z) (2*M2(r,z) - sqrt(4*M2(r,z)^2 - 4*(n2(r,z)-omega0^2)*(xi2(r,z)-omega0^2)))/(2*(n2(r,z)-omega0^2));

    cgr = @(r,z,l,m) (m/(omega0*(l^2+m^2)^2))*(n2(r,z)*l*m + M2(r,z)*(l^2-m^2) - l*m*xi2(r,z));
    cgz = @(r,z,l,m) (-l/(omega0*(l^2+m^2)^2))*(n2(r,z)*l*m + M2(r,z)*(l^2-m^2) - l*m*xi2(r,z));
    
    domegar = @(r,z,l,m) (1/(2*omega0*(l^2+m^2)))*(dn2r(r,z)*l^2 - 2*m*l*dM2r(r,z) + m^2*dxi2r(r,z));
    domegaz = @(r,z,l,m) (1/(2*omega0*(l^2+m^2)))*(dn2z(r,z)*l^2 - 2*m*l*dM2z(r,z) + m^2*dxi2z(r,z));

    %starting values
    y0 = [rstart; zstart; m0*(ar1(rstart,zstart)*(ch==1) + ar2(rstart,zstart)*(ch==2)); m0]; tspan = [0 t];
    r_ray(1) = y0(1); z_ray(1) = y0(2); l_ray(1) = y0(3); m_ray(1) = y0(4);tmat(1)  = 0; tend = tmat(end);

    %ODE's to solve
    raytraceODEs = @(tau,y) [cgr(y(1),y(2),y(3),y(4));cgz(y(1),y(2),y(3),y(4));-domegar(y(1),y(2),y(3),y(4));-domegaz(y(1),y(2),y(3),y(4))];
    
    %if the well boundary is encounter
    options = odeset('RelTol',err,'AbsTol',err,'Events',@lookout);

    while tend < t
        %ode solving phase
        [tau, y, ~, ~,ie] = ode45(raytraceODEs,tspan,y0,options);
    
        %append values
        nl = length(tau);
        tmat = [tmat;tau(2:end-1)];
        r_ray = [r_ray; y(2:end-1,1)];
        z_ray = [z_ray; y(2:end-1,2)];
        l_ray = [l_ray; y(2:end-1,3)];
        m_ray = [m_ray; y(2:end-1,4)];
        tend = tau(end);

        %Calling back?
        if ((sum(ie==1) + sum(ie==3) + sum(ie==4))>0) && (tau(end)<t)
            disp('touched boundaries');
            break;
        elseif ((sum(ie==2))>0) && (tau(end)<t)
            disp(['reflecting at:: ',num2str(tau(end))]);
            tspan = [tau(end-1) t];
            if ch==1
                ch = 2;
            else
                ch = 1;
            end
            y0 = [y(end-1,1),y(end-1,2), y(end-1,4)*(ar1(y(end-1,1),y(end-1,2))*(ch==1) + ar2(y(end-1,1),y(end-1,2))*(ch==2)),y(end-1,4)];
        elseif tau(end)>= t
            disp('completed');
            break;
        end
    end
    for i=1:1:length(tmat)
        cgz_ray(i,1) = cgz(r_ray(i),z_ray(i),l_ray(i),m_ray(i));
        cgr_ray(i,1) = cgr(r_ray(i),z_ray(i),l_ray(i),m_ray(i));
        omegamat(i,1) = omega(r_ray(i),z_ray(i),l_ray(i),m_ray(i));
    end
    
    function [position, isterminal, direction] = lookout(tau,y)
        position = [y(1)-Lr; y(1); y(2)-Lz; y(2);];
        isterminal = [1; 1; 1; 1;];
        direction = [-1; -1; -1; -1;];
    end
end
