clear
%f - Coriolis frequency, z0 - initial location, H - mixed layer depth
%nb - N (Background Brunt-Vaisala frequency), h - e-folding depth
%vm - maximum velocity, rm - radius at maximum velocity, b - radius of the
%shield divided by rm, a - power of the radial function in CLM vortex

%bri = 15; bcu = 0.1;
f = 1*10^(-4); z0 = 0; b=4; a = 2;
H = 700; h = H/2; 
vm = -1; rm = 25000;
Lr = 7*rm; Lz = -2*H;
nb = 50*f; bcu = 2*vm/(f*rm);

%(rstart, zstart) - starting point of ray, t - maximum time to run
%dt - timestep, Lr - radial extent of domain, Lz - vertical extent of the
%domain, omega0 - frequency of the ray, m0 - initial vertical wavenumber
%ch - characteristic to consider (1 or 2). 
niwdays = 10^5; ip = 2*pi/f; t = niwdays*ip;
rstart = 1.3*rm; zstart = -0.03*H;
omega0 = 0.9999*f; lambdaz = 100; m0 = 2*pi/lambdaz;


%calculating the minimum frequency contours
%n2 = Brunt-Vaisala frequency, m2 - baroclinic term, xi2 - modified Coriolis
%frequency, ro - Rossby no., ri - Richardson no., cu - Curvature no., 
% wmin1 and wmin2 - two roots of minimum frequency
[n2, m2, xi2, rhomat, ri, ro, cu, wmin1, wmin2]  = minfreq_CLM(a, z0, H, f, nb,vm, rm, Lr, Lz);


if sum(sum(wmin1/f<1)) > 0
    disp(['wmin1 can also have trapping? at r_m = ', num2str(rm),', vm =',num2str(vm)]);
end

figobj=figure(1);
hold on
figobj.Position = [10 10 700 500];
rmat = (1/100:1/100:Lr/rm)*rm; zmat = (-1/100:-1/100:Lz/H)*H;
[X,Y] = meshgrid(rmat/rm,zmat/H);
wmmat = wmin2; wmmat(imag(wmmat)~=0) = NaN;
wmmatmin = (min(real(wmmat(:,1)/f))); disp(wmmatmin);
[C0,d0] = contour(X,Y,rhomat',15,'Color',[.7 .7 .7]);
%hold on
%d0.LineWidth = 1;
[C,d] = contour(X,Y,real(wmmat)'/f,[wmmatmin, 0.9999, 0.9999],'k');
d.LineWidth = 2;
clabel(C,d); fontsize(25,"points");
hold on
[C,d] = contour(X,Y,real(wmin2.^2)',[0,  0],'b'); d.LineWidth = 1;
xlabel("$r/r_m$"); ylabel("$z/H$");
title(['$Cu_b =$ ',num2str(bcu)]);
xlim([0 5]); ylim([-2 0]);
box on

%%

%main location. (r_ray,z_ray) - location of the ray, (l_ray, m_ray) -
%radial and vertical wavenumbers of the ray,(cgr_ray, cgz_ray) - radial and 
%vertical group velocities along the ray, chmat - characteristic along the ray

ch = 1;
[r_ray1,z_ray1,l_ray1,m_ray1,tmat1,omegamat1,cgr_ray1,cgz_ray1,chmat1] = raytrace_CLM_ode23(rstart,zstart,t,Lr,Lz,omega0,m0,ch,a,vm,rm,z0,h,f,nb);
cg1 = sqrt(real(cgr_ray1).^2 + real(cgz_ray1).^2);
ch = 2;
[r_ray2,z_ray2,l_ray2,m_ray2,tmat2,omegamat2,cgr_ray2,cgz_ray2,chmat2] = raytrace_CLM_ode23(rstart,zstart,t,Lr,Lz,omega0,m0,ch,a,vm,rm,z0,h,f,nb);
cg2 = sqrt(real(cgr_ray2).^2 + real(cgz_ray2).^2);
hold on
plot(r_ray1((imag(r_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8) & (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/rm,z_ray1((imag(z_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8)& (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/H,'r', 'LineWidth',1);
hold on
plot(r_ray2((imag(r_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/rm,z_ray2((imag(z_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/H,'g', 'LineWidth',1);
xlim([0 5]); ylim([-2 0]);

% figobj=figure(2);
% figobj.Position = [10 10 900 900];
% subplot(3,1,1);
% nmat = (n2')*H^2/(f^2*rm^2);
% [C,d] = contourf(X,Y,nmat,100); set(d,'LineColor','None'); 
% colorbar; clim([min(min(nmat)),max(max(nmat))]);
% hold on
% plot(r_ray1((imag(r_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8) & (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/rm,z_ray1((imag(z_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8)& (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/H,'r', 'LineWidth',2);
% hold on
% plot(r_ray2((imag(r_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/rm,z_ray2((imag(z_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/H,'g', 'LineWidth',2);
% hold on
% [C2,d2] = contour(X,Y,real(wmmat)'/f,[omega0/f, omega0/f],'k'); clabel(C2,d2); 
% fontsize(14,"points");
% [C3,d3] = contour(X,Y,(wmin2.^2)',[0,  0],'b'); d3.LineWidth = 2;
% xlabel("$r/r_m$"); ylabel("$z/H$");
% title(['$N^2H^2/f^2L^2$ ($N_b/f =$ ',num2str(nb/f), '; $v_m =$ ',num2str(vm),' m/s; $r_m = $ ',num2str(rm),' m;)']);
% xlim([0 5]); ylim([-2 0]);
% 
% subplot(3,1,2);
% [C,d] = contourf(X,Y,(xi2')/f^2,100); set(d,'LineColor','None'); 
% colorbar; clim([min(min(xi2/f^2)),max(max(xi2/f^2))]);
% hold on
% plot(r_ray1((imag(r_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8) & (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/rm,z_ray1((imag(z_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8)& (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/H,'r', 'LineWidth',2);
% hold on
% plot(r_ray2((imag(r_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/rm,z_ray2((imag(z_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/H,'g', 'LineWidth',2);
% hold on
% [C2,d2] = contour(X,Y,real(wmmat)'/f,[omega0/f, omega0/f],'k'); clabel(C2,d2); 
% fontsize(14,"points");
% [C3,d3] = contour(X,Y,(wmin2.^2)',[0,  0],'b'); d3.LineWidth = 2;
% xlabel("$r/r_m$"); ylabel("$z/H$");
% title(['$\chi^2/f^2$ ($N_b/f =$ ',num2str(nb/f), '; $v_m =$ ',num2str(vm),' m/s; $r_m = $ ',num2str(rm),' m;)']);
% xlim([0 5]); ylim([-2 0]);
% 
% subplot(3,1,3);
% mmat = m2'*H/(f^2*rm);
% [C,d] = contourf(X,Y,mmat,100); set(d,'LineColor','None'); 
% colorbar; clim([min(min(mmat)),max(max(mmat))]);
% hold on
% plot(r_ray1((imag(r_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8) & (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/rm,z_ray1((imag(z_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8)& (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/H,'r', 'LineWidth',2);
% hold on
% plot(r_ray2((imag(r_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/rm,z_ray2((imag(z_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/H,'g', 'LineWidth',2);
% hold on
% [C2,d2] = contour(X,Y,real(wmmat)'/f,[omega0/f, omega0/f],'k'); clabel(C2,d2); 
% fontsize(14,"points");
% [C3,d3] = contour(X,Y,(wmin2.^2)',[0,  0],'b'); d3.LineWidth = 2;
% xlabel("$r/r_m$"); ylabel("$z/H$");
% title(['$M^2H/f^2L$ ($N_b/f =$ ',num2str(nb/f), '; $v_m =$ ',num2str(vm),' m/s; $r_m = $ ',num2str(rm),' m;)']);
% xlim([0 5]); ylim([-2 0]);


% figobj=figure(2);
% figobj.Position = [10 10 900 900];
% qd = (1+ro)-(1+cu).*ri.^(-1);
% ld = (1+cu); pd = qd.*ld;
% subplot(3,1,1);
% [C,d] = contourf(X,Y,ld',100); set(d,'LineColor','None'); 
% colorbar; clim([min(min(ld)),max(max(ld))]);
% hold on
% plot(r_ray1((imag(r_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8) & (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/rm,z_ray1((imag(z_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8)& (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/H,'r', 'LineWidth',2);
% hold on
% plot(r_ray2((imag(r_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/rm,z_ray2((imag(z_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/H,'g', 'LineWidth',2);
% hold on
% [C2,d2] = contour(X,Y,real(wmmat)'/f,[omega0/f, omega0/f],'k'); clabel(C2,d2); 
% fontsize(14,"points");
% [C3,d3] = contour(X,Y,(wmin2.^2)',[0,  0],'b'); d3.LineWidth = 2;
% xlabel("$r/r_m$"); ylabel("$z/H$");
% title(['$L = (1+Cu)$ ($N_b/f =$ ',num2str(nb/f), '; $v_m =$ ',num2str(vm),' m/s; $r_m = $ ',num2str(rm),' m;)']);
% xlim([0 5]); ylim([-2 0]);
% 
% subplot(3,1,2);
% %minri = min(min(abs(ri)));
% %riplot = ri; riplot(abs(ri)>(minri+10)) = NaN;
% [C,d] = contourf(X,Y,qd',100); set(d,'LineColor','None'); 
% colorbar; clim([min(min(qd)),max(max(qd))]);
% hold on
% plot(r_ray1((imag(r_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8) & (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/rm,z_ray1((imag(z_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8)& (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/H,'r', 'LineWidth',2);
% hold on
% plot(r_ray2((imag(r_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/rm,z_ray2((imag(z_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/H,'g', 'LineWidth',2);
% hold on
% [C2,d2] = contour(X,Y,real(wmmat)'/f,[omega0/f, omega0/f],'k'); clabel(C2,d2); 
% fontsize(14,"points");
% [C3,d3] = contour(X,Y,(wmin2.^2)',[0,  0],'b'); d3.LineWidth = 2;
% xlabel("$r/r_m$"); ylabel("$z/H$");
% title(['$q = (1+Ro)-(1+Cu)Ri^{-1}$ ($N_b/f =$ ',num2str(nb/f), '; $v_m =$ ',num2str(vm),' m/s; $r_m = $ ',num2str(rm),' m;)']);% Create textbox
% xlim([0 5]); ylim([-2 0]);
% % annotation(figure(3),'textbox',...
% %     [0.849080133555927 0.633220910623946 0.100836393989983 0.0311973018549739],...
% %     'String','$\log_{10}(Ri)$',...
% %     'FontSize',18,...
% %     'FitBoxToText','off',...
% %     'EdgeColor','none');

% subplot(2,1,1);
% pd2 = pd; pd2((pd>omega0/f)) = NaN;
% [C,d] = contourf(X,Y,pd',100); set(d,'LineColor','None'); 
% colorbar; clim([min(min(pd)),1]);
% hold on
% plot(r_ray1((imag(r_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8) & (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/rm,z_ray1((imag(z_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8)& (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/H,'r', 'LineWidth',2);
% hold on
% plot(r_ray2((imag(r_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/rm,z_ray2((imag(z_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/H,'g', 'LineWidth',2);
% hold on
% [C2,d2] = contour(X,Y,real(wmmat)'/f,[omega0/f, omega0/f],'k'); clabel(C2,d2); 
% fontsize(14,"points");
% [C3,d3] = contour(X,Y,(wmin2.^2)',[0,  0],'b'); d3.LineWidth = 2;
% xlabel("$r/r_m$"); ylabel("$z/H$");
% title('$\Phi = Lq$');
% xlim([0 5]); ylim([-2 0]);
% fontsize(25,'points');
% 
% subplot(2,1,2);
% [C,d] = contourf(X,Y,(xi2')/f^2,100); set(d,'LineColor','None'); 
% colorbar; clim([min(min(xi2/f^2)),1]);
% hold on
% plot(r_ray1((imag(r_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8) & (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/rm,z_ray1((imag(z_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8)& (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/H,'r', 'LineWidth',2);
% hold on
% plot(r_ray2((imag(r_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/rm,z_ray2((imag(z_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/H,'g', 'LineWidth',2);
% hold on
% [C2,d2] = contour(X,Y,real(wmmat)'/f,[omega0/f, omega0/f],'k'); clabel(C2,d2); 
% fontsize(14,"points");
% [C3,d3] = contour(X,Y,(wmin2.^2)',[0,  0],'b'); d3.LineWidth = 2;
% xlabel("$r/r_m$"); ylabel("$z/H$");
% title('$\chi^2/f^2$ ');
% xlim([0 5]); ylim([-2 0]);
% fontsize(25,'points');

%figobj = figure(2);
%[C,d] = contourf(X,Y,sqrt(n2'/f^2),100); set(d,'LineColor','None'); 
%colorbar; clim([min(min(sqrt(n2/f^2))),max(max(sqrt(n2/f^2)))]);
%hold on
%plot(r_ray1((imag(r_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8) & (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/rm,z_ray1((imag(z_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8)& (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/H,'r', 'LineWidth',2);
%hold on
%plot(r_ray2((imag(r_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/rm,z_ray2((imag(z_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/H,'g', 'LineWidth',2);
%hold on
%[C2,d2] = contour(X,Y,real(wmmat)'/f,[omega0/f, omega0/f],'k'); clabel(C2,d2); 
%fontsize(14,"points");
%[C3,d3] = contour(X,Y,real(wmin2.^2)',[0,  0],'b'); d3.LineWidth = 2;
%xlabel("$r/r_m$"); ylabel("$z/H$");
%title('$N^2/f^2$ ');
%xlim([0 5]); ylim([-2 0]);
%fontsize(25,'points');

%figobj = figure(3);
%figobj.Position = [10 10 900 900];
%[C,d] = contourf(X,Y,(xi2')/f^2,100); set(d,'LineColor','None'); 
%colorbar; clim([min(min(xi2/f^2)),1]);
%hold on
%plot(r_ray1((imag(r_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8) & (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/rm,z_ray1((imag(z_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8)& (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/H,'r', 'LineWidth',2);
%hold on
%plot(r_ray2((imag(r_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/rm,z_ray2((imag(z_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/H,'g', 'LineWidth',2);
%hold on
%[C2,d2] = contour(X,Y,real(wmmat)'/f,[omega0/f, omega0/f],'k'); clabel(C2,d2); 
%fontsize(14,"points");
%[C3,d3] = contour(X,Y,real(wmin2.^2)',[0,  0],'b'); d3.LineWidth = 2;
%xlabel("$r/r_m$"); ylabel("$z/H$");
%title('$\chi^2/f^2$ ');
%xlim([0 5]); ylim([-2 0]);
%fontsize(25,'points');

% figobj = figure(3);
% %figobj.Position = [10 10 900 900];
% wmmat2 = wmmat; wmmat2(wmmat/f>=1) = NaN;
% [C,d] = contourf(X,Y,real(wmmat2')/f,100); set(d,'LineColor','None'); 
% colorbar; clim([min(min(real(wmmat2)/f)),1]);
% hold on
% plot(r_ray1((imag(r_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8) & (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/rm,z_ray1((imag(z_ray1)==0) & (sqrt(cgr_ray1.^2+cgz_ray1.^2)>1e-8)& (imag(cgr_ray1)==0) & (imag(cgz_ray1)==0))/H,'r', 'LineWidth',2);
% hold on
% plot(r_ray2((imag(r_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/rm,z_ray2((imag(z_ray2)==0) & (sqrt(cgr_ray2.^2+cgz_ray2.^2)>1e-8)& (imag(cgr_ray2)==0) & (imag(cgz_ray2)==0))/H,'g', 'LineWidth',2);
% hold on
% [C2,d2] = contour(X,Y,real(wmmat)'/f,[omega0/f, omega0/f],'k'); clabel(C2,d2); 
% fontsize(14,"points");
% [C3,d3] = contour(X,Y,(wmin2.^2)',[0,  0],'b'); d3.LineWidth = 2;
% xlabel("$r/r_m$"); ylabel("$z/H$");
% title('$\omega_min/f$ ');
% xlim([0 5]); ylim([-2 0]);
% fontsize(25,'points');
