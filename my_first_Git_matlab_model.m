clear
close all

Ix=100;
Lx=1;
x=linspace(0,Lx,Ix); x=x(1:end-1);
u = 0.5+1*sin(2*pi*x/Lx);
delta = 0.0001;
dx=x(2)-x(1);
T = 2;
dt=0.7*min(0.5*dx^2/delta,0.5*dx/max(u));
t=0;

figure;
plot(x,u);
while t<T
   t=t+dt;
   dt=0.7*min(0.5*dx^2/delta,0.5*dx/max(u));
   u=timestep(u,dt,dx,delta);
   plot(x,u);
   xlim([0 1]*Lx);
   ylim([-1 1]*1.5);
   title(['t=',num2str(t,'%2.2f')]);
   drawnow;
end

function u=timestep(u,dt,dx,delta)
    ui   = u;
    uip1 = circshift(u,-1);
    uim1 = circshift(u,1 );
    ri = (ui-uim1)./(uip1-ui);
    rim1 = circshift(ri,1);
    phi_i = fluxlim_vanleer(ri);
    phi_im1 = fluxlim_minmod(rim1);
    uLimh = uim1+0.5*phi_im1.*(ui-uim1);
    uRimh = ui-0.5*phi_i.*(uip1-ui);
    Simh = 0.5*(uLimh+uRimh);
    
    uimh_star = zeros(size(uLimh));    
    uimh_star(uLimh>uRimh & Simh>0)   = uLimh(uLimh>uRimh & Simh>0);
    uimh_star(uLimh>uRimh & Simh<0)   = uRimh(uLimh>uRimh & Simh<0);
    uimh_star(uLimh==uRimh)           = uLimh(uLimh==uRimh);
    uimh_star(uLimh<uRimh & uLimh>0)  = uLimh(uLimh<uRimh & uLimh>0);
    uimh_star(uLimh<uRimh & uRimh<0)  = uRimh(uLimh<uRimh & uRimh<0);  
    uimh_star(uLimh<0 & uRimh>0)      = 0;
   
    Fimh = uimh_star.^2;
    Fiph = circshift(Fimh,-1);

    dudt = -(Fiph-Fimh)/dx + delta/dx^2*(-2*u+circshift(u,1)+circshift(u,-1));  % scheme using a flux limiter;
%     dudt = -(ui.^2-uim1.^2)/dx + 0*delta/dx^2*(-2*u+circshift(u,1)+circshift(u,-1));
    u=u+dt*dudt;
end

function phi = fluxlim_minmod(r)
    phi=max(0,min(1,r));
    phi(isinf(r))=1;
end

function phi = fluxlim_vanleer(r)
    phi=(r+abs(r))./(1+abs(r));
    phi(isinf(r))=2;
end
