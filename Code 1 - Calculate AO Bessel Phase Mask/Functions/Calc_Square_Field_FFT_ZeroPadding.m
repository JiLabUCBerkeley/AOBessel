function [x_Fieldout,y_Fieldout,FieldOut]=Calc_Square_Field_FFT_ZeroPadding(FieldIn,x_FieldIn,y_FieldIn,f1,lambda,FFT2D_padding)
%lambda, mm;
%f1, mm;
nVector=(1:1:FFT2D_padding)-ceil(FFT2D_padding/2);
ImgCal_temp=zeros(FFT2D_padding,FFT2D_padding);
ImgCal_temp((end-length(x_FieldIn))/2+1:(end+length(x_FieldIn))/2,(end-length(y_FieldIn))/2+1:(end+length(y_FieldIn))/2)=FieldIn;
step_x=mean(diff(x_FieldIn));
step_y=mean(diff(y_FieldIn));

x_Fieldout=(lambda*1e3*f1*1e3)/(FFT2D_padding*step_x).*nVector;%um fx=(lambda*f)/(N.deltax)*n
y_Fieldout=(lambda*1e3*f1*1e3)/(FFT2D_padding*step_y).*nVector;%um
% x_Fieldout(floor(end/2)+1)=[];
% y_Fieldout(floor(end/2)+1)=[];
FieldOut=fftshift(fft2(fftshift(ImgCal_temp)));
end
