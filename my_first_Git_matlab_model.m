clear
close all

Ix=100;
Lx=1;
x=linspace(0,Lx,Ix); x=x(1:end-1);
u = 1+0.1*sin(2*pi*x/Lx);
delta = 0.02;
dx=x(2)-x(1);
T = 1;
dt=min(0.25*dx^2/delta,0.5*dx/max(u));
It=round(T/dt);

figure;
plot(x,u);
for it=1:It
   u=timestep(u,dt,dx,delta);
   plot(x,u);
   xlim([0 1]*Lx);
   ylim([0 1]*1.5);
   drawnow;
end

function u=timestep(u,dt,dx,delta)
    dudt = (u.^2-circshift(u,-1).^2)/dx + delta/dx^2*(-2*u+circshift(u,1)+circshift(u,-1));
    u=u+dt*dudt;
end

