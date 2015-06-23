

re = -0.5 * ones(1,10);  im= linspace(-0.2,0.2,10);
scatter(re,im);
xlim([-1,1]);
ylim([-1,1]);

m = re + 1i * im;
plot(angle(m));

gradient = [1,0,-1];

% Generate Test Data
Row = 64; Col = 64;
img = zeros(64,64)
theta = linspace(0, 6 * pi, Row);
for r = 1:Row
    for c = 1:Col
        img(r,c) = angle(exp(i* theta(r)));
    end
end

figure();
imshow(img,[]);

% Sobel 
Gx = [1, 0, -1; 2 0 -2; 1 0 -1];
Gy = [1 2 1; 0 0 0; -1 -2 -1];
Gx = conv2(Gx,img);
Gy = conv2(Gy,img);
G = sqrt(Gx.*Gx + Gy.*Gy);

figure();
imshow(G,[]);

