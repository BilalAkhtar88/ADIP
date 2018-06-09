function tileBlur()

clc;
clear all;

I = imread('cameraman.tif');

if numel(size(I)) >= 3
    I = rgb2gray(I);
end

G = double(I);
M = size(G,1);
N = size(G,2);
alpha = 10;
sigma = 10;
% h = fspecial('gaussian',60,sigma);
% size(G)
Iconv = imgaussfilt(G,sigma);

for i = 1:alpha
    Aprime(i,:) = Iconv(M - alpha + i,:);
    Aprime(M + alpha + i, :) = Iconv(i,:);
    Bprime(:,i) = Iconv(:,N - alpha + i);
    Bprime(:,N + alpha + i, :) = Iconv(:,i);
end

mask = Aprime < 1;
mask = imfill(mask,'holes');
Aprime = regionfill(Aprime, mask);

mask = Bprime < 1;
mask = imfill(mask,'holes');
Bprime = regionfill(Bprime, mask);


% xA_co = [1 N N 1];
% yA_ro = [alpha alpha M+alpha+1 M+alpha+1];
% xB_co = [alpha alpha+N alpha+N alpha];
% yB_ro = [1 1 M M];


for i = 1:alpha
    Cprime(i,:) = Bprime(M-alpha+i,:);
    Cprime(M+alpha+i,:) = Bprime(i,:);
end
for i = 1:alpha
    Cprime(:,i) = Aprime(:,N-alpha+i);
    Cprime(:,N+alpha+i) = Aprime(:,i);
end

mask = Cprime < 1;
mask = imfill(mask,'holes');
Cprime = regionfill(Cprime, mask);

% xC_co = [alpha alpha+N alpha+N alpha];
% yC_ro = [alpha alpha alpha+M alpha+M];
% 
% Cprime = regionfill(Cprime, xC_co, yC_ro);

A = Aprime(alpha+1:alpha+M,1:N);
B = Bprime(1:M,alpha+1:alpha+N);
C = Cprime(alpha+1:alpha+M,alpha+1:alpha+N);

bigIm(1:M,1:N) = C;
bigIm(1:M,N+1:2*N) = A;
bigIm(1:M,2*N+1:3*N) = C;
bigIm(M+1:2*M,1:N) = B;
bigIm(M+1:2*M,N+1:2*N) = G;
bigIm(M+1:2*M,2*N+1:3*N) = B;
bigIm(2*M+1:3*M,1:N) = C;
bigIm(2*M+1:3*M,N+1:2*N) = A;
bigIm(2*M+1:3*M,2*N+1:3*N) = C;



figure(2)
subplot(3,1,1)
imshow(Aprime,[])
subplot(3,1,2)
imshow(Bprime,[])
subplot(3,1,3)
imshow(Cprime,[])

figure(3)
imshow(bigIm,[])

end