function zernikeCoeff = ZernikeDecomposition(wavefront, centerPos, pupilSize,ZernikeModeN)

% wavefront = wavefront - mean(wavefront(:));
wfInPupil = wavefront(centerPos(2)-pupilSize/2+1:centerPos(2)+pupilSize/2,...
                    centerPos(1)-pupilSize/2+1: centerPos(1)+pupilSize/2);
ImSize = size(wfInPupil, 1);
[xx, yy] = meshgrid(1:ImSize);
xx = xx - pupilSize/2;
yy = yy - pupilSize/2;
% wfInPupil(xx.^2+yy.^2 > pupilSize^2/4) = 0;%%Turnate the corners
% figure,imagesc(wfInPupil);axis image;

% ZernikeModeN = 55;
zernikeCoeff = zeros(1, ZernikeModeN);
for i = 1:ZernikeModeN
    [Mode, ~] = zernike_fun(i, ImSize);                                % Mode: 51x51 wavefront for each Zernike mode
    zernikeCoeff(i) = sum(wfInPupil(:) .* Mode(:)) / sum(Mode(:).^2);   % Zernike_coef*mode=phase --> Zer_coef=phase/mode 
end

end