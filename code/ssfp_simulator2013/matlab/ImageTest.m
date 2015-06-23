clear; clf;

% Create a Gaussian
[x, y] = meshgrid(-4:0.05:4);
s1 = exp(-pi * (x.^2 + y.^2));
imshow(s1, []);

% Create a 2D Sinc
s2 = sin(-pi*x).*sin(-pi*y)./(pi^2*x.*y);
imshow(s2, []);

% Create a 2D Sinc without lines
[x, y] = meshgrid(-4:0.1:4 + 10*eps);
s2 = sin(-pi*x).*sin(-pi*y)./(pi^2*x.*y);
imshow(s2, []);

% Exponential Signal
kx0 = 0.25;
ky0 = 0.5;
s3 = exp(1i*2*pi*(kx0*x + ky0*y));
subplot(1,2,1);
imshow(real(s3),[]);
subplot(1,2,2);
imshow(imag(s3),[]);

% Create a Rect
clf;
s4 = double(abs(x)<0.5 & abs(y)<0.5);
imshow(s4,[]);

% Load Face Image
f1 = imread('shepp256.png');
f1 = double(f1);
imshow(f1, []);

% Load Tissue Image 64x64
% [img,map] = imread('Phantom64Tissues.png');
% figure, image(img), axis equal, colormap(map)
f2 = imread('Phantom64Tissues.png');
imshow(f2, []);

% black = [0 0 0];
% white = [255 255 255];
% red = [237 28 36];

% FFT Practice
s = f1;
S = fftshift(fft2(s));
os = ifft2(S);
subplot(1,3,1);
imshow(abs(s),[]);       % Orignal Signal
subplot(1,3,2);
imshow(log(abs(S)),[]);  % Log Mag of FFT of Signal
subplot(1,3,3); 
% imshow(abs(os),[]);    % iFFT(FFT of Signal)
imshow(angle(S),[]);     % Angle of FFT of Signal
