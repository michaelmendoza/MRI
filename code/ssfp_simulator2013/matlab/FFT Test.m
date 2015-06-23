x = [1 2 3 4 5 5 4 3 2 1];
N = length(x);
n = 0:N-1;
X = zeros(N,1);
for k = 0:N-1
% 	X(k+1) = sum(x .* exp(-1i*2*pi*k*n/N));
    X(k+1) = sum(x .* exp(1i*2*pi*k*n/N));
end
fftshift(X)

x = [1 2 3 4 5 5 4 3 2 1];
N = length(x);
n = 0:N-1;
X = zeros(N,1);
for k = -N/2:N/2 - 1 
% 	X(k+1) = sum(x .* exp(-1i*2*pi*k*n/N));
    X(k+N/2+1) = sum(x .* exp(1i*2*pi*k*n/N));
end
X