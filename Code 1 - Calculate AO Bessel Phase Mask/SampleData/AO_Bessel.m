%% AO-Bessel.m
% The AO-Bessel.m was designed to generate an aberration-corrected Bessel phase mask for
% high-efficency Bessel AO correction at the objective focal plane.
%%
clear all;
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));
%% System parameters:
global SLM;
SLM.pixelNumber=512;
SLM.pitch=15;
SLM.dimention=SLM.pitch*SLM.pixelNumber;

global lambda;
lambda=940e-6;
global Magnification_SLM2toMask;
Magnification_SLM2toMask=1/(150/30*350/750);
global f1;
f1=200;
global beamD_SLM1;
beamD_SLM1=2;

%% Mask
global Mask;
Mask.innerDiameter =1.015;
Mask.outerDiameter=1.2;
global S;
S=22.5;
PhasePattern_SLM1=double(AnnularGratingOnSLMGenerator(S,SLM.pixelNumber,SLM.pixelNumber))./256*2*pi;%rad
step_SLM1=8; 

R_SLM1=0:step_SLM1:SLM.dimention/2;
Amplitude_SLM1=1*exp(-R_SLM1.^2/(beamD_SLM1*1000/2)^2); 
Phase_SLM1=interp1(SLM.pitch*(0:1:SLM.pixelNumber/2),PhasePattern_SLM1(SLM.pixelNumber/2, SLM.pixelNumber/2:end),R_SLM1);%rad
Field_SML1_0=Amplitude_SLM1.*exp(1i*Phase_SLM1);
 
step_Mask=3;
D_mask=4;
R_Mask=0:step_Mask:D_mask*1000/2;
Field_Mask_0=Fourier_CircularLens(Field_SML1_0,lambda*1e3,R_SLM1,f1,R_Mask);

Field_Mask2D_Amp=Creat3DPSF(R_Mask,abs(Field_Mask_0));

Field_Mask2D_Phase=Creat3DPSF(R_Mask,angle(Field_Mask_0));
Field_MaskSection_Amp=[rot90(Field_Mask2D_Amp,2)  rot90(Field_Mask2D_Amp,1); rot90(Field_Mask2D_Amp,3),Field_Mask2D_Amp];
Field_MaskSection_Amp(:,end/2)=[];
Field_MaskSection_Amp(end/2,:)=[];

Field_MaskSection_Phase=[rot90(Field_Mask2D_Phase,2)  rot90(Field_Mask2D_Phase,1); rot90(Field_Mask2D_Phase,3),Field_Mask2D_Phase];
Field_MaskSection_Phase(:,end/2)=[];
Field_MaskSection_Phase(end/2,:)=[];


x_RingMask=[-flip(R_Mask,2),R_Mask(2:end)];
y_RingMask=[-flip(R_Mask,2),R_Mask(2:end)];

Amplitude_RingMask=Field_MaskSection_Amp;
Phase_temp=Creat3DPSF(R_Mask,angle(Field_Mask_0));
Phase_RingMask=[rot90(Phase_temp,2)  rot90(Phase_temp,1); rot90(Phase_temp,3),Phase_temp];
Phase_RingMask(:,end/2)=[];
Phase_RingMask(end/2,:)=[];
Field_RingMask=Amplitude_RingMask.*exp(1i.*Phase_RingMask);
Field_RingMask(isnan(Field_RingMask))=0;

titleSize=12;
fig12=figure(1);
ax1=subplot(2,2,1);
imagesc(x_RingMask,y_RingMask,Field_MaskSection_Amp);
h_title=title('Ideal-Bessel Amplitude at Mask');h_title.FontSize=titleSize;
xlabel('x (um)');
ylabel('y (um)');
colormap(ax1,jet)
colorbar;

ax2=subplot(2,2,2);
imagesc(x_RingMask,y_RingMask,Field_MaskSection_Phase);
h_title=title('Ideal-Bessel Phase at Mask');h_title.FontSize=titleSize;
xlabel('x (um)');
ylabel('y (um)');
colormap(ax1,jet)
h1 = colorbar;
set(get(h1,'title'),'string','Phase (rad)');

ax3=subplot(2,2,3);
plot(x_RingMask,Field_MaskSection_Amp(ceil(end/2),:),'k','linewidth',1.5);
h_title=title('Ideal-Bessel Amplitude Profile');h_title.FontSize=titleSize;
xlabel('x (um)');
ylabel('Amplitude (a.u.)');
vline([-Mask.outerDiameter*1e3/2,-Mask.innerDiameter*1e3/2,Mask.outerDiameter*1e3/2,Mask.innerDiameter*1e3/2]);


ax4=subplot(2,2,4);
plot(x_RingMask,Field_MaskSection_Phase(ceil(end/2),:),'k','linewidth',2);
h_title=title('Ideal-Bessel Phase Profile');h_title.FontSize=titleSize;
xlabel('x (um)');
ylabel('Optical field phasee (rad)');
vline([-Mask.outerDiameter*1e3/2,-Mask.innerDiameter*1e3/2,Mask.outerDiameter*1e3/2,Mask.innerDiameter*1e3/2]);


%% Load AO pattern
[AOfile,AOFilePath] = uigetfile('*.tif');
AOPhase=imread([AOFilePath AOfile],'tif');
AOPhase=double(AOPhase/255*2*pi);
[AOImSizeX,AOImSizeY]=size(AOPhase);

AOPhaseMap=zeros(512,512);
AOPhaseMap(512/2-AOImSizeX/2+1:512/2+AOImSizeX/2,512/2-AOImSizeY/2+1:512/2+AOImSizeY/2)=AOPhase;
AOPhaseMap=imresize(AOPhaseMap,[512 512],'bilinear');

R_SLM2=3;
x_SLM2=(0:R_SLM2:SLM.dimention)-SLM.dimention/2;
y_SLM2=(0:R_SLM2:SLM.dimention)-SLM.dimention/2;

AOPhaseMap2=imresize(AOPhaseMap,[length(x_SLM2) length(y_SLM2)],'cubic');
Amplitude_SLM2=ones([length(x_SLM2) length(y_SLM2)]);
Field_SLM2=Amplitude_SLM2.*exp(1i.*AOPhaseMap2);

ConjugateTimes=3;
if(mod(ConjugateTimes,2))
Field_MaskAO=rot90(Field_SLM2,2);
else
Field_MaskAO=Field_SLM2;    
end

x_MaskAO=x_SLM2.*Magnification_SLM2toMask;%um
y_MaskAO=y_SLM2.*Magnification_SLM2toMask;%um 

[xTepIn,yTepIn]= meshgrid(x_RingMask,y_RingMask);
[xTepOut,yTepOut]= meshgrid(x_MaskAO,y_MaskAO);
Field_RingMask_ReSampling=interp2(xTepIn,yTepIn,Field_RingMask,xTepOut,yTepOut);

R_MaskAO=x_MaskAO(end/2:end);
step_MaskAO=mean(diff(x_MaskAO));
[~,~,Mask2D_MaskAO_Amp,f_MaskAO]=Generate2DAnnularApodizingMask(R_MaskAO,Mask.innerDiameter,Mask.outerDiameter);

Amplitude_Masked_RingAO=Mask2D_MaskAO_Amp;
phase_Masked_RingAO=rot90(AOPhaseMap2,2).*Mask2D_MaskAO_Amp;
Field_Masked_RingAO=Amplitude_Masked_RingAO.*exp(1i.*phase_Masked_RingAO);

FFT2D_padding=8192*4;
[x_SLM1Back,y_SLM1Back,Field_SLM1Back_tem]=Calc_Square_Field_FFT_ZeroPadding(Field_Masked_RingAO,x_MaskAO,y_MaskAO,f1,lambda,FFT2D_padding);
x_SLM1Back(or(x_SLM1Back<-SLM.dimention/2,x_SLM1Back>SLM.dimention/2))=[];
y_SLM1Back(or(y_SLM1Back<-SLM.dimention/2,y_SLM1Back>SLM.dimention/2))=[];
Field_SLM1Back=Field_SLM1Back_tem((end-length(x_SLM1Back))/2+1:(end+length(x_SLM1Back))/2,(end-length(y_SLM1Back))/2+1:(end+length(y_SLM1Back))/2);
Amplitude_SLM1Back=abs(Field_SLM1Back);
Phase_SLM1Back=imgaussfilt(angle(Field_SLM1Back),0.8);
clear Field_SLM1Back_tem;

PhasePattern_new=imresize(Phase_SLM1Back,[SLM.pixelNumber,SLM.pixelNumber],'bilinear');
AmplitudePattern_new=imresize(Amplitude_SLM1Back,[SLM.pixelNumber,SLM.pixelNumber],'bilinear');
X_SLM1=-SLM.dimention/2:step_SLM1:SLM.dimention/2;
Y_SLM1=-SLM.dimention/2:step_SLM1:SLM.dimention/2;
[x_SLM1_new,y_SLM1_new,Field_new]=GenerateVirtualSLMGausianFiled(PhasePattern_new,beamD_SLM1,X_SLM1,Y_SLM1);

FFT2D_padding=4096*4;
[X_MaskConfirm,Y_MaskConfirm,OutField_tem]=Calc_Square_Field_FFT_ZeroPadding(Field_new,x_SLM1_new,y_SLM1_new,f1,lambda,FFT2D_padding);
X_MaskConfirm(or(X_MaskConfirm<-SLM.dimention/2,X_MaskConfirm>SLM.dimention/2))=[];
Y_MaskConfirm(or(Y_MaskConfirm<-SLM.dimention/2,Y_MaskConfirm>SLM.dimention/2))=[];
OutField=OutField_tem((end-length(X_MaskConfirm))/2+1:(end+length(X_MaskConfirm))/2,(end-length(Y_MaskConfirm))/2+1:(end+length(Y_MaskConfirm))/2);
OutField=rot90(OutField,2);
clear OutField_tem;

Amplitude_MaskConfirm=abs(OutField);
Phase_MaskConfirm=imgaussfilt(angle(OutField),0.8);

Hfig17=figure(2);
ax1=subplot(2,2,1);
imagesc(X_MaskConfirm,Y_MaskConfirm,Amplitude_MaskConfirm);
xlim([-2000, 2000]);
ylim([-2000, 2000]);
title('AO-Bessel Amplitude at Mask');

ax2=subplot(2,2,2);
imagesc(X_MaskConfirm,Y_MaskConfirm,Phase_MaskConfirm);
xlim([-2000, 2000]);
ylim([-2000, 2000]);
caxis([-pi,pi]);
colormap(ax2, jet);
h1 = colorbar;
set(get(h1,'title'),'string','Phase (rad)');
title('AO-Bessel Phase at Mask');

ax3=subplot(2,2,3);
plot(X_MaskConfirm,Amplitude_MaskConfirm(ceil(end/2),:),'k','linewidth',1.5);
title('AO-Bessel Amplitude Profile');
xlabel('x (um)');
xlim([-2000, 2000]);
ylabel('Amplitude (a.u.)');
vline([-Mask.outerDiameter*1e3/2,-Mask.innerDiameter*1e3/2,Mask.outerDiameter*1e3/2,Mask.innerDiameter*1e3/2]);

ax4=subplot(2,2,4);
plot(X_MaskConfirm,Phase_MaskConfirm(ceil(end/2),:),'g','linewidth',1);
hold on;
plot(x_MaskAO,phase_Masked_RingAO(ceil(end/2),:),'b','linewidth',1.5);

title('AO-Bessel Phase Profile');
xlabel('x (um)');
xlim([-1500, 1500]);
ylabel('Phase (rad)');
vline([-Mask.outerDiameter*1e3/2,-Mask.innerDiameter*1e3/2,Mask.outerDiameter*1e3/2,Mask.innerDiameter*1e3/2]);

Hfig18=figure(3);
axs1=subplot(1,2,1);
imagesc(x_MaskAO,y_MaskAO,phase_Masked_RingAO);
xlim([-1500, 1500]);
ylim([-1500, 1500]);
xlabel('um');
ylabel('um');
colormap(axs1,jet);
caxis([0 pi]);
h1 = colorbar;
set(get(h1,'title'),'string','Phase (rad)');
title('Measured Phase within the Ring');
hold on;
plot([-fliplr(R_MaskAO),R_MaskAO],[fliplr(f_MaskAO'),f_MaskAO'].*R_MaskAO(end),'r','linewidth',2);

R_MaskConfirm=X_MaskConfirm(floor(end/2)+1:end);
[~,~,MaskOutConfirm,MaskOutConfirm_Profile]=Generate2DAnnularApodizingMask(R_MaskConfirm,Mask.innerDiameter,Mask.outerDiameter);
axs2=subplot(1,2,2);
imagesc(X_MaskConfirm,Y_MaskConfirm,medfilt2(Phase_MaskConfirm,[2,2]).*MaskOutConfirm);
xlim([-1500, 1500]);
ylim([-1500, 1500]);
xlabel('um');
ylabel('um');
colormap(axs2,jet);
caxis([0 pi]);
h1 = colorbar;
set(get(h1,'title'),'string','Phase (rad)');
title('Generated Phase within the Ring');
hold on;
plot([-fliplr(R_MaskConfirm),R_MaskConfirm],[fliplr(MaskOutConfirm_Profile'),MaskOutConfirm_Profile'].*R_MaskConfirm(end),'r','linewidth',2);

axs19=figure(4);
imagesc(1:1:SLM.pixelNumber,1:1:SLM.pixelNumber,PhasePattern_new);
xlabel('SLM1 Pixel index (1~512)');
ylabel('SLM1 Pixel index (1~512)');
xlim([1,SLM.pixelNumber]);
ylim([1,SLM.pixelNumber]);
h1 = colorbar;
set(get(h1,'title'),'string','Phase (rad)');
colormap(axs19, jet);
caxis([-pi pi]);

title('AO-Bessel Phase Pattern');

%% Save AO_Phase Mask SLM1
PhasePattern_newSave=(PhasePattern_new/pi+1)./2*255;
filename='AO_Bessel_Focal';
imwrite(uint8(PhasePattern_newSave),[filename '.bmp']);