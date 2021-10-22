function [PSF] = Calc_Annular_Field_Integrals_V2(x, y, z,field, field_x, field_y, wavelength,refind,f0)
%CALC_ANNULAR_FIELD_INTEGRALS  Calculates PSF for non-rotational-symmetric phase pattern
%on SLM
%based on "Principeles of nano-optics 3.51"
%f0:mm, wavelength:um, xyz:um, field_x/y:mm

f=f0*refind*1e3;  %um focal length in media
k=2*pi/wavelength*refind; %um^-1
[field_xx,field_yy]=meshgrid(field_x,field_y);  %mm
field_rr=sqrt(field_xx.^2+field_yy.^2)+eps;
field_th=atan(field_rr./f*1e3);
field_phi=angle(field_xx+1.j*field_yy);
dx=(max(field_x)-min(field_x))/(length(field_x)-1);
dy=(max(field_y)-min(field_y))/(length(field_y)-1);

[field_inf_x,field_inf_y,field_inf_z]=Calc_field_inf(field,field_th, field_phi,refind);

Ex=zeros(length(x),length(y),length(z));
Ey=zeros(length(x),length(y),length(z));
Ez=zeros(length(x),length(y),length(z));
x_len=length(x);
y_len=length(y);
z_len=length(z);

for a=1:x_len
    for b=1:y_len
        for c=1:z_len
            rho=sqrt(x(a).^2+y(b).^2)+eps;
            psi=angle(x(a)+1.j*y(b));
            integrand_x=field_inf_x.*exp(1.j.*k.*z(c).*cos(field_th)+1.j*k.*rho.*sin(field_th).*cos(field_phi-psi)).*sin(field_th);
            integrand_y=field_inf_y.*exp(1.j.*k.*z(c).*cos(field_th)+1.j*k.*rho.*sin(field_th).*cos(field_phi-psi)).*sin(field_th);
            integrand_z=field_inf_z.*exp(1.j.*k.*z(c).*cos(field_th)+1.j*k.*rho.*sin(field_th).*cos(field_phi-psi)).*sin(field_th);
            integrand_x=integrand_x*(1.j*k*f*exp(-1.j*k*f))/2/pi;
            integrand_y=integrand_y*(1.j*k*f*exp(-1.j*k*f))/2/pi;
            integrand_z=integrand_z*(1.j*k*f*exp(-1.j*k*f))/2/pi;
            %polar coordinate to xy coordinate. use sine condition th=r/f
            integrand_x=integrand_x*f*1e-3./(field_rr.*(field_rr.^2+(f*1e-3).^2));
            integrand_y=integrand_y*f*1e-3./(field_rr.*(field_rr.^2+(f*1e-3).^2));
            integrand_z=integrand_z*f*1e-3./(field_rr.*(field_rr.^2+(f*1e-3).^2));
%             Ex(a,b,c)=sum(integrand_x,'all')*dx*dy;
%             Ey(a,b,c)=sum(integrand_y,'all')*dx*dy;
%             Ez(a,b,c)=sum(integrand_z,'all')*dx*dy;
            Ex(a,b,c)=sum(sum(integrand_x,'omitnan'),'omitnan')*dx*dy;
            Ey(a,b,c)=sum(sum(integrand_y,'omitnan'),'omitnan')*dx*dy;
            Ez(a,b,c)=sum(sum(integrand_z,'omitnan'),'omitnan')*dx*dy;
        end
    end
end

%%

PSF=abs(Ex).^2+abs(Ey).^2+abs(Ez).^2;
PSF=PSF.^2;
end