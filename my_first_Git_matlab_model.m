clear
close all

fluxlimtype = 'vanleer'; % 'vanleer'/'minmod'/'superbee';
Ix=100;
Lx=1;
x=linspace(0,Lx,Ix); x=x(1:end-1);
u = 0.5+1*sin(2*pi*x/Lx);
delta = 0.001;
dx=x(2)-x(1);
T = 2;
dt=0.7*min(0.5*dx^2/delta,0.5*dx/max(u));
t=0;

figure;
plot(x,u);
while t<T
   t=t+dt;
   dt=0.3*min(0.5*dx^2/delta,0.5*dx/max(u));
   u=timestep(u,dt,dx,delta,fluxlimtype);
   plot(x,u);
   xlim([0 1]*Lx);
   ylim([-1 1]*1.5);
   title(['t=',num2str(t,'%2.2f')]);
   drawnow;
end

function u=timestep(u,dt,dx,delta,fluxlimtype)
    ui   = u;
    uip1 = circshift(u,-1);
    uim1 = circshift(u,1 );
    ri = (ui-uim1)./(uip1-ui);
    rim1 = circshift(ri,1);
    switch fluxlimtype
        case 'vanleer'
            phi_i   = fluxlim_vanleer(ri);
            phi_im1 = fluxlim_vanleer(rim1); 
        case 'minmod'
            phi_i   = fluxlim_minmod(ri);
            phi_im1 = fluxlim_minmod(rim1);
        case 'superbee'
            phi_i   = fluxlim_superbee(ri);
            phi_im1 = fluxlim_superbee(rim1);
        otherwise
            error(['Limiter type: ',fluxlimtype,' does not exist!']);
    end
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

function phi = fluxlim_superbee(r)
    phi=max([0,min(2*r,1),min(r,2)]);
    phi(isinf(r) & sign(r)==1)=2;
    phi(isinf(r) & sign(r)==-1)=0;
end
