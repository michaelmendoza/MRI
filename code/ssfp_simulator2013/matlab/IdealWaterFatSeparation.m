% Water Fat Separation
% Michael Mendoza
 
clear;
%%
% 2 Point Dixon Reconstruction
if(0)
 
% m=load('.\Stephen\IMS_606_607.mat');
% m1=m.im1; m2=m.im2;
% m=load('.\Stephen\IMS_632_633.mat');
% m1=squeeze(m.im1(100,:,:)); m2=squeeze(m.im2(100,:,:));
 
m1amp=dicomread('.\Breast_Anonymous\gm_gre te 4.7 - 10\IM-0001-0100.dcm'); 
m1phase=dicomread('.\Breast_Anonymous\gm_gre te 4.7 - 11\IM-0001-0100.dcm');
m1amp=cast(m1amp,'single'); m1phase=cast(m1phase,'single');
m1phase=(m1phase/4095)*2*pi;
m1=m1amp.*exp(i*m1phase);
 
m2amp=dicomread('.\Breast_Anonymous\gm_gre te 5.75 - 8\IM-0001-0100.dcm');
m2phase=dicomread('.\Breast_Anonymous\gm_gre te 5.75 - 9\IM-0001-0100.dcm');
m2amp=cast(m2amp,'single'); m2phase=cast(m2phase,'single');
m2phase=(m2phase/4095)*2*pi;
m2=m2amp.*exp(i*m2phase);
 
mw=0.5*(m1+m2);
mf=0.5*(m1-m2);
 
figure(1);
%imshow(abs(mw)',[0 0.07])
imshow(abs(mw)',[0 460])
figure(2);
%imshow(abs(mf)',[0 0.07])
imshow(abs(mf)',[0 660])
 
end
%%
% 3 Point Dixon Reconstruction
if(0)
 
m1amp=dicomread('.\Breast_Anonymous\gm_gre te 4.7 - 10\IM-0001-0100.dcm'); 
m1phase=dicomread('.\Breast_Anonymous\gm_gre te 4.7 - 11\IM-0001-0100.dcm');
m1amp=cast(m1amp,'single'); m1phase=cast(m1phase,'single');
m1phase=(m1phase/4095)*2*pi;
m1=m1amp.*exp(i*m1phase);
 
m2amp=dicomread('.\Breast_Anonymous\gm_gre te 5.75 - 8\IM-0001-0100.dcm');
m2phase=dicomread('.\Breast_Anonymous\gm_gre te 5.75 - 9\IM-0001-0100.dcm');
m2amp=cast(m2amp,'single'); m2phase=cast(m2phase,'single');
m2phase=(m2phase/4095)*2*pi;
m2=m2amp.*exp(i*m2phase);
 
m3amp=dicomread('.\Breast_Anonymous\gm_gre te 6.8 - 6\IM-0001-0100.dcm');
m3phase=dicomread('.\Breast_Anonymous\gm_gre te 6.8 - 7\IM-0001-0100.dcm');
m3amp=cast(m3amp,'single'); m3phase=cast(m3phase,'single');
m3phase=(m3phase/4095)*2*pi;
m3=m3amp.*exp(i*m3phase);
 
phi=angle(conj(m1).*m3)/2;
 
mw3=0.5*(m1+m2.*exp(-i*phi));
mf3=0.5*(m1-m2.*exp(-i*phi));
 
figure(3);
imshow(abs(mw3)',[0 460]);
figure(4);
imshow(abs(mf3)',[0 660]);
 
end
 
%% Comparsion
if(0)
    figure(5);
    imshow(abs(mw3)'-abs(mw)',[0 460/10]);
    figure(6);
    imshow(abs(mf3)'-abs(mf)',[0 660/10]);
end
 
%%
% IDEAL
 
if(1)
Threshold=0.5;
dF=[0 430];
%dF=[0 -220];
 
% Load Data
load idealexample.mat;
t=[-0.0004 0.0011 0.0027];
%t =  [-0.000378787878788 0.001136363636364 0.002651515151515];
%t=[0.0004 0.0011 0.0027];
m1=IM1; m2=IM2; m3=IM3; 
N = 3; M = 2;  %Number of Images/ Number of Species
im(1,:,:)=m1;
im(2,:,:)=m2;
im(3,:,:)=m3;
% load IMS_606_607.mat;
% 
% imshow(IMS_606_607,[])
 
% Load S vector
J=length(m1(:,1));
K=length(m1(1,:));
Sdata=zeros(2*N,J,K);
S_dhat=zeros(6,1);
gR=zeros(1,N);
gI=zeros(1,N);
B=zeros(2*N,5);
dpR=zeros(1,M);
dpI=zeros(1,M);
%p =zeros(6,1);
pimage1=zeros(J,K);
pimage2=zeros(J,K);
for n=1:N
    for j=1:J
        for k=1:K
            Sdata(n,j,k)  =real( im(n,j,k) );
            Sdata(N+n,j,k)=imag( im(n,j,k));
        end
    end
end
 
% Fill c and d
for j=1:M % Species
   for k=1:N % Images
      c(j,k)=cos(2*pi*dF(j)*t(k));
      d(j,k)=sin(2*pi*dF(j)*t(k)); 
   end
end
 
% Fill A 
j=2; % Using special case for A in Appendix C
A=zeros(2*N,4);
for k=1:N
    A(k,:)    =[1, 0, c(j,k), -d(j,k)];
    A(N+k,:)  =[0, 1, d(j,k), c(j,k) ];
end
 
% Loop through for each pixel
for px=1:K
%for px=110:140
    for py=1:J     
        
        %Load Data
        S_hat=Sdata(:,px,py);
        
        % Estimate singal of each species
        p = (A.'*A)\(A.'*S_hat);
        
        % Set Flags and Counters
        wflag=1;
        count=0;
        
        % Check if pixel is zero
        if(Sdata(1,px,py)==0 && Sdata(2,px,py)==0)            
           wflag=0; 
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Loop through while psi error is > threshold
        dPsi = 2*Threshold;
        psi0 = 0;
        while(dPsi>Threshold && wflag==1)
                     
            % Fill B 
            j=2; % Using special case for B in Appendix C
            for k=1:N
                gR(k)=2*pi*t(k)*(-p(2)-p(3)*d(j,k)-p(4)*c(j,k));
                gI(k)=2*pi*t(k)*(p(1)+p(3)*c(j,k)-p(4)*d(j,k));
                B(k,:)    =[gR(k), 1, 0, c(j,k), -d(j,k)];
                B(N+k,:)  =[gI(k), 0, 1, d(j,k), c(j,k) ];
            end
 
            % Estimate of error
            % Calc S_dhat
            for n=1:N
                temp0=0;
                temp1=0;
                for j=1:M
                    idR=2*(j-1)+1;
                    idI=2*(j-1)+2;
                    temp0=temp0+(p(idR)*c(j,n)-p(idI)*d(j,n));
                    temp1=temp1+(p(idR)*d(j,n)+p(idI)*c(j,n));
                end
                S_dhat(n)  =S_hat(n)-temp0;
                S_dhat(N+n)=S_hat(N+n)-temp1;
            end
 
    
            [px,py]
            % Calculate Error Estimates
            y=(B.'*B)\(B.'*S_dhat);
 
            dPsi=y(1);
            for j=1:M
                dpR(j)=y(2*j);
                dpI(j)=y(2*j+1);
            end
 
            % Recalculate Psi
            psi0=psi0+dPsi; 
    
            % Recalculate S_hat
            for n=1:N
                Stemp=(Sdata(n,px,py)+i*Sdata(n+N,px,py))*exp(-i*2*pi*psi0*t(n));
                S_hat(n)  =real(Stemp);
                S_hat(n+N)=imag(Stemp);
            end
    
            count=count+1
            psivec(count) = psi0;
            dPsivec(count)= dPsi;
        end
        
        % Save Fat and Water Consentrations
        if(wflag==1)
            pimage1(px,py)=sqrt(p(1)*p(1)+p(2)*p(2));
            pimage2(px,py)=sqrt(p(3)*p(3)+p(4)*p(4));
        else
            pimage1(px,py)=0;
            pimage2(px,py)=0;
        end 
    end
end
 
figure(1);
imshow(pimage1(:,:),[]);
figure(2);
imshow(pimage2(:,:),[]);
figure(3);
imshow(abs(m1),[]);
end
