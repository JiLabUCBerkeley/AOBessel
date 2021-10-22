function [Ex,Ey,Ez]=Calc_field_inf(field,th, phi,refind)
%Calculate electric field at infinite distance. based on "Principeles of nano-optics 3.51"
Ex=field/2.*((1+cos(th))-(1-cos(th)).*cos(2*phi)).*sqrt(1/refind).*sqrt(cos(th));
Ey=field/2.*(cos(th)-1).*sin(2*phi).*sqrt(1/refind).*sqrt(cos(th));
Ez=field/2.*(-2).*cos(phi).*sin(th).*sqrt(1/refind).*sqrt(cos(th));
end