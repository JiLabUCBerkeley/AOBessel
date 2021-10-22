function [xMask,yMask,MaskOut,MaskOut_Profile]=Generate2DAnnularApodizingMask(R_Mask,Mask_innerDiameter,Mask_outerDiameter)
% R_Mask,um, [0:step:end] the radical coordinates, e.g., R_Mask=x_MaskAO(end/2:end);
% Mask_innerDiameter,mm
% Mask_outerDiameter, mm
% xMask, um
% yMask, um
% MaskOut, double 1 or 0;
% MaskOut_Profile, Cross sectional profile MaskOut_Profile



step_Mask=mean(diff(R_Mask));
f_Mask=zeros(length(R_Mask),1);
f_Mask(round(Mask_innerDiameter*1e3/step_Mask/2):round(Mask_outerDiameter*1e3/step_Mask/2))=1;

Mask2D_Mask=Creat3DPSF(R_Mask,f_Mask);%phase in unit pi
Mask2D_Mask(isnan(Mask2D_Mask))=0;
Mask2D_MaskAO_Amp=[rot90(Mask2D_Mask,2)  rot90(Mask2D_Mask,1); rot90(Mask2D_Mask,3),Mask2D_Mask];
Mask2D_MaskAO_Amp(:,end/2)=[];
Mask2D_MaskAO_Amp(end/2,:)=[];

xMask=[-fliplr(R_Mask),R_Mask(1:end)];
yMask=[-fliplr(R_Mask),R_Mask(1:end)];
MaskOut=Mask2D_MaskAO_Amp;
MaskOut_Profile=f_Mask;
end