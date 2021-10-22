%% Creat a 2D Gaussian
function mat = Generate2DGaussianSurface(mat,BeamSize,center)
%%
% mat --- dimention of the squire field, in um
% mat --- dimention of beamsize, in mm

%%
gsize= mat;
for r=1:gsize(1)
    for c=1:gsize(2)
        mat(r,c) = gaussC(r,c, BeamSize,center);
    end
end

end

function val = gaussC(x, y, BeamSize,center)
 xc = center(1);
 yc = center(2);
exponent = ((x-xc).^2 + (y-yc).^2)./(BeamSize*1000/2)^2;
val       = (exp(-exponent));    
end