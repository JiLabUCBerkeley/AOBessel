% Fourier transfrom for an ideal lens. Input field is placed at the back
% focal plane and the output field is at the focal plane.
function field_output=Fourier_CircularLens(field_input,wavelength,r_input,f,r_output)
% f-mm, r_input, output in um
f=f*1000; % f in um
p=length(r_output);  
field_output=zeros(size(r_output));
for ii=1:p
    a=field_input.*besselj(0,2*pi*r_input*r_output(ii)/wavelength/f)*2*pi.*r_input;
    field_output(ii)=trapz(r_input,a);
end
end
