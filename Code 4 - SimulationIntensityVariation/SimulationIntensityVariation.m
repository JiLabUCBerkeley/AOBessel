clear all;
%Simulation mode:
flag_Bessel=1;
flag_Gaussian=0;
%Error mode:
flag_addPhaseError=0;
%Correction mode:
flag_focalAO=1;
flag_pupilAO=0;

Abb_Z=200*1e3;

%% SLM parameters2
SLM.pitch=15; %um
SLM.pixelNumber=512;
SLM.size=SLM.pitch*SLM.pixelNumber; %um
upper_bound=255;  %gray level range
lower_bound=0;
delta_degree=2*pi/(upper_bound-lower_bound);
laser_sig=3.262; %mm
refind=1.33;
lambda=940*1e-3; %um
obj.mag=25;
obj.NA=1.05;
obj.f=180/obj.mag; %mm in air
obj.D=2*obj.f*obj.NA; %mm
mag=2; %magnification
k=2*pi/lambda; %um^-1
SlitX1=5107;
SlitX2=5580;
SlitX10=5207;
SlitX20=5480;
currentFolder=pwd;
Resultpath=[currentFolder '\FocalAO\'];
load([currentFolder '\field_int.mat']);


%% field at back pupil plane
phase1=imread([currentFolder '\AstiAO_Cal.tif']);
phase1=pi/255*double(phase1)-pi;

phase=ones(size(phase1));
phase=2*pi/255*double(phase)-pi;
[pixel_num,~]=size(phase);


%coordinate
field_x=(-pixel_num/2:pixel_num/2-1)*SLM.pitch*2;
field_y=(-pixel_num/2:pixel_num/2-1)*SLM.pitch*2; 
[field_xx,field_yy]=meshgrid(field_x,field_y);
field_rr=sqrt(field_xx.^2+field_yy.^2);

%% Uniform intensity
intensity=zeros(pixel_num);
if(flag_Bessel)
intensity(field_rr>SlitX1&field_rr<SlitX2)=1;
else
intensity(:)=1;%Gaussian
end
field=intensity.*exp(1.j*phase);

%% Propogate field and add phase error other than the pupil
zeroPadding=4;
intensity_new=zeros(4*size(intensity));
phase_new=zeros(4*size(phase));
phase_error=zeros(4*size(phase));

intensity_new((end-length(intensity))/2+1:(end+length(intensity))/2,(end-length(intensity))/2+1:(end+length(intensity))/2)=intensity;
phase_new((end-length(phase))/2+1:(end+length(phase))/2,(end-length(phase))/2+1:(end+length(phase))/2)=phase;
phase_error((end-length(phase1))/2+1:(end+length(phase1))/2,(end-length(phase1))/2+1:(end+length(phase1))/2)=phase1;

N_upSample=4096*2;
phase_error=imresize(phase_error,[N_upSample N_upSample]);
field_upSample_Intensity=imresize(intensity_new,[N_upSample N_upSample]);
field_upSample_Phase=imresize(phase_new,[N_upSample N_upSample]);
field_upSample=field_upSample_Intensity.*exp(1.j*field_upSample_Phase);
field_temp=propIR(field_upSample,field_rr(end/2,end)*2,lambda,Abb_Z);
%%Add phase error
if(flag_addPhaseError)
phase_Z=angle(field_temp)+phase_error;
else
phase_Z=angle(field_temp);
end
intensity_Z=abs(field_temp);
field_temp_Z=intensity_Z.*exp(1.j*phase_Z);

field_temp2=propIR(field_temp_Z,field_rr(end/2,end)*2,lambda,-Abb_Z);
field_back=field_temp2(length(field_temp2)*3/8+1:length(field_temp2)*5/8,length(field_temp2)*3/8+1:length(field_temp2)*5/8);

figure(1);
imagesc(AB);
colormap(jet);
title('With focal AO');

figure(2);
subplot(2,2,1);
imagesc(field_x,field_y,abs(field));
title('Amplitude without aberration');
subplot(2,2,2);
imagesc(field_x,field_y,angle(field));
title('Phase without aberration');
x=-1023:1:1024;
y=-1023:1:1024;
[XX,YY]=meshgrid(x,y);
rr=sqrt(XX.^2+YY.^2);
rr_mask=ones(size(rr));
rr_mask(rr>760)=0;
rr_mask(rr<670)=0;
subplot(2,2,3);
int_back=abs(field_back);
int_back(~rr_mask)=0;
imagesc(field_x,field_y,int_back)
title('Amplitude with aberration');
subplot(2,2,4);
phase_back=angle(field_back);
phase_back(~rr_mask)=0;
imagesc(field_x,field_y,phase_back)
title('Phase with aberration');
%% Construc the final field
int_tem=imresize(abs(field_back),[512,512]);
int_tem1=zeros([512,512]);
int_tem2=zeros([512,512]);
int_tem3=zeros([512,512]);
int_tem1(field_rr>SlitX1&field_rr<SlitX2)=int_tem(field_rr>SlitX1&field_rr<SlitX2);
int_tem2(field_rr>SlitX1&field_rr<SlitX2)=AB(field_rr>SlitX1&field_rr<SlitX2);
if(flag_addPhaseError)
    if(flag_focalAO)
AB1=AB*sum(sum(sum(int_tem1^2))/sum(sum(int_tem2^2))).*int_tem2;
    else
AB1=int_tem2;        
    end
else
AB1=int_tem1;    
end
int_tem3(field_rr>SlitX10&field_rr<SlitX20)=AB1(field_rr>SlitX10&field_rr<SlitX20);
%%

if(flag_Gaussian)
 field=int_tem.*exp(1.j*1);   
else
 field=int_tem3.*exp(1.j*1);   
end

field_x2=field_x;
field_y2=field_y;
field2=field;
if(flag_Bessel)
x=-2:0.05:2; %Bessel
y=-2:0.05:2;
z=-40:5:40;
else
x=-1:0.05:1; %Gaussian
y=-1:0.05:1;
z=-2:0.2:2;
end

for a=1:length(z)
    PSF=Calc_Annular_Field_Integrals_V2(x, y, z(a),field2, field_x2*1e-3, field_y2*1e-3, lambda,refind,obj.f);
    PSF=uint16(squeeze(PSF)*1e-9);
    imwrite(PSF,[Resultpath,num2str(z(a)),'.tif']);
end

img=imread([Resultpath,num2str(z(1)),'.tif']);
imwrite(img, [Resultpath,'Stack_0.05xy_5z.tif'], 'tif', 'WriteMode', 'overwrite','Compression', 'none');

for a=2:length(z)
    img=imread([Resultpath,num2str(z(a)),'.tif']);
    imwrite(img, [Resultpath,'Stack_0.05xy_5z.tif'], 'tif', 'WriteMode', 'append','Compression', 'none');
end

















