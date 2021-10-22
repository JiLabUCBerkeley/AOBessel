function [Z Mask] = zernike_fun(ModeN,ZSize)
% ----------
% r: radial distance;   
% n, Rm: integer indeices, m<=n, n>=0
% R: R(n,m) radial polynomials
% Sequetial indices: OSA and ANSI single-index Zernike polynomials (j=n(n-2)+m/2)
% ModeN: j+1, starting from 1
% ----------
%% ModeN >= 1;

x = linspace(-1, 1, ZSize);
[x, y] = meshgrid(x, x);
r = sqrt(x.^2+y.^2);
theta = acos(x./(r+eps)); % eps for 0/0 case
theta = theta.*(y>=0) - theta.*(y<0);
Mask = (r <= 1);%Circular
% Mask=ones(size(r));%sqaure
% hint: j-1=(n(n+2)+m)/2 >=(n^2+n)/2 --> 8(j-1)+1>=4n^2+4n+1
        %  --> n<=(sqrt(8(j-1)+1)-1)/2
        %  --> m = (j-1)*2-n(n+2) - n(n-1)
n = floor((sqrt(8*(ModeN-1)+1)-1) /2);
m = (ModeN-1)*2-n*(n+2);
% n1 = ceil((sqrt(1+8*ModeN)-1)/2); % Raph's code
% m1 = ModeN-n1*(n1-1)/2;
% n = n1-1;
% m = -n+(m1-1)*2;
Rm = abs(m);
R = 0;

for k = 0:(n-Rm)/2
    R = R + (-1)^k *factorial(n-k) / ...
        (factorial(k)*factorial((n+Rm)/2-k)*factorial((n-Rm)/2-k)) ...
        * r.^(n-2*k);
end

if m >=0
    Z = R.*cos(Rm*theta);
else
    Z = R.*sin(Rm*theta);
end

Z = Z.*Mask;

%% rms normalization
avr = mean(Z(:));
% rms = sqrt(sum((Z(:)-avr).^2) / numel(Z));
rms = sqrt(sum((Z(:)-avr).^2) / sum(Mask(:)));
Z = Z / rms;

% fprintf('ModeN=%d, n=%d, m=%d, rms =%f\n', ModeN, n, m, rms); % Show indices