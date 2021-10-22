function phasePattern=AnnularGratingOnSLMGenerator(S,pixelNumberX,pixelNumberY)
% patternOnSLMGenerator is to generate binary pattern used for SLM;
% S is the period of binary pattern on SLM;
% pixelNumberX is the dimension of SLM on x axis;
% pixelNumberY is the dimension of SLM on y axis;
%
%---------example------------------
% S=22.4;
% pixelNumberX=512;
% pixelNumberY=512;
% phasePattern=patternOnSLMGenerator(S,pixelNumberX,pixelNumberY);   

if nargin==0
    S=22.5;
    pixelNumberX=512;
    pixelNumberY=512;
end
m=1:pixelNumberY;
n=1:pixelNumberX;
[nn,mm]=meshgrid(n,m);
yy=1*(mm-round(pixelNumberY/2));
xx=1*(nn-round(pixelNumberX/2));
rr=yy.^2+xx.^2;rr=sqrt(rr);
beta=1*pi/4;
phasePattern_tmp=cos(2*pi*rr/S-beta)+eps;
phasePattern=zeros(pixelNumberY,pixelNumberX);
phasePattern(phasePattern_tmp<0)=0;
phasePattern(phasePattern_tmp>0)=1;
phasePattern=uint8(phasePattern*256/2);
end
