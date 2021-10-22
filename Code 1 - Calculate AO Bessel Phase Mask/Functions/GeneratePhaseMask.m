function PhaseMask=GeneratePhaseMask(type,rBO,rPhase,xPixShift,yPixShift)
%% Parameters
if nargin==0
type='circular';
rBO=0:0.4:7600;%%radius at the back pupil plan,um
rPhase=mod(linspace(0, 80, length(rBO)),pi)'./pi*2-1;%phase at back pupil plane
xPixShift=0;
yPixShift=0;
end
mp=2;%Magnification from SLM to backpupil plane
% rPhase=-rPhase;% Flip the phase value
rBO_ResizeToSLM=rBO./mp;% Mask coordinate at SLM   
pixelSize=15;%um
rSLM=0:pixelSize:pixelSize*512/2-1;
rPhaseSLM=interp1(rBO_ResizeToSLM,rPhase,rSLM,'next');


% 
% figure(1)
% plot(rBO_ResizeToSLM,rPhase,'*-r');
% hold on;
% plot(rSLM,rPhaseSLM,'o-k');

%% Map and rezie circular mask to SLM
if(type=='circular')
PhaseValue=zeros(512,512);
xC=256+xPixShift;
yC=256+yPixShift;
    for i1=1:1:512
       for i2=1:1:512
           rIndex=round(sqrt(((i1-xC)^2+(i2-yC)^2)));
           if(rIndex<=256&&rIndex~=0)
       PhaseValue(i1,i2)=rPhaseSLM(rIndex); 
       PhaseValue(i1,i2)=(PhaseValue(i1,i2)+1)/2*255;
           else
       PhaseValue(i1,i2)=0;
           end
       end
    end
end
PhaseMask=uint8(PhaseValue);
%% Save image mask
figure(1);
imshow(PhaseMask);
% SaveMaskPath = uigetdir;
% imwrite(PhaseMask,[SaveMaskPath '\Mask.png']);

end