function [X_SLM1,Y_SLM1,OutField]=GenerateVirtualSLMGausianFiled(PhasePattern_SLM1,beamD_SLM1,X_SLM1,Y_SLM1)

SLM_pitch=15;%um
OutField=zeros(length(X_SLM1),length(X_SLM1));

Phase_SLM1=imresize(PhasePattern_SLM1,[length(X_SLM1),length(Y_SLM1)]);
beamD_SLM1_Nor=beamD_SLM1/SLM_pitch;
Amplitude_SLM1=Generate2DGaussianSurface([512, 512],beamD_SLM1_Nor,[512/2, 512/2]);
Amplitude_SLM1=imresize(Amplitude_SLM1,[length(X_SLM1),length(Y_SLM1)]);
OutField=Amplitude_SLM1.*exp(1i*Phase_SLM1);

end