function tileBlur()

clc;
clear all;
%% Reading Input Image

I = imread('cameraman.tif');
if numel(size(I)) >= 3
    I = rgb2gray(I);
end

figure(1)
subplot(2,3,1)
imshow(I,[])
title('Original Image')

M = size(I,1);
N = size(I,2);

pSig1 = max(max(I(1:floor(0.1*M),:)));
pSig2 = max(max(I(size(I,1)-floor(0.1*M):size(I,1),:)));
pSig3 = max(max(I(:,size(I,2)-floor(0.1*N):size(I,2))));
pSig4 = max(max(I(:,1:floor(0.1*N))));
pSig = max([pSig1 pSig2 pSig3 pSig4]);

numIt = 20;

%% Adding blur to the input Image

sigmaBlur = 1;
len = 10;

PSF = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 1 1 1 1 1 1 0 0 0 0 0 0 0;
    0 0 0 0 0 0 1 1 0 0 0 1 1 1 1 1 0 0 0 0;
    0 0 0 0 0 1 1 1 0 0 0 1 1 1 1 0 0 0 0 0;
    0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 0 0 0 0 0;
    0 0 0 0 0 0 1 1 1 1 0 0 1 1 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 0;
    0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 0 0;
    0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
    ];
PSF = PSF*0.01;
% PSF = fspecial('gaussian', len, sigmaBlur); %Uniform gaussian blur for experiment
G = imfilter(I, PSF, 'conv');
figure(1)
subplot(2,3,2)
imshow(G,[])
title('Result of convolution with non-uniform PSF')

%% Deconvolving the blurred Image without solving for Boundary Value Problem

J1 = deconvlucy(G,PSF,numIt);
figure(2)
imshow(PSF,[])
figure(1)
subplot(2,3,4)
imshow(J1,[])

numElem = numel(J1(1:floor(0.1*M),:)) + numel(J1(:,1:floor(0.1*N))) + numel(J1(:,size(G,2)-floor(0.1*N):size(G,2))) + numel(J1(size(G,1)-floor(0.1*M):size(G,1),:)); 
noiseSqr = sqrt(sum(sum((J1(1:floor(0.1*M),:)-I(1:floor(0.1*M),:)).^2)));
noiseSqr = noiseSqr + sqrt(sum(sum((J1(:,1:floor(0.1*N))-I(:,1:floor(0.1*N))).^2)));
noiseSqr = noiseSqr + sqrt(sum(sum((J1(:,size(G,2)-floor(0.1*N):size(G,2))-I(:,size(G,2)-floor(0.1*N):size(G,2))).^2)));
noiseSqr = noiseSqr + sqrt(sum(sum((J1(size(G,1)-floor(0.1*M):size(G,1),:)-I(size(G,1)-floor(0.1*M):size(G,1),:)).^2)));
noiseSqr = noiseSqr / numElem;

pSNR = db(double(pSig)/noiseSqr)
hT1 = title('Restored Image without filling boundaries');
T1Pos = round(get(hT1,'Position')); %// Get the position
hT1_2 = text(T1Pos(1),T1Pos(2) + size(J1,2)+10,strcat('PSNR =', {' '}, num2str(pSNR,'%.2f'), 'dB'),'HorizontalAlignment','center'); %// Place the text


%% Deconvolving after using edgetaper(.) command

G = G(len+1:size(G,1)-len,len+1:size(G,2)-len);
% G = G(1:size(G,1)-len,1:size(G,2)-len); %Cropping G matrix of the boundary values that need to be filled; 
%VVV IMPORTANT!!!!!!!!!!!!!!!!!!!!! 
%This cropping information is missing from the paper

M = size(G,1);
N = size(G,2);
tic
% timerVal = tic
G = edgetaper(G,PSF);
J2 = deconvlucy(G,PSF,numIt);
toc
% elapsedTime = toc
% toc(timerVal)
% elapsedTime = toc(timerVal)

figure(1)
subplot(2,3,5)
imshow(J2,[])
hT2 = title('Restored Image after using edgetaper(.)')
T2Pos = round(get(hT2,'Position')); %// Get the position

I = I(len+1:size(I,1)-len,len+1:size(I,2)-len);

pSig1 = max(max(I(1:floor(0.1*M),:)));
pSig2 = max(max(I(size(I,1)-floor(0.1*M):size(I,1),:)));
pSig3 = max(max(I(:,size(I,2)-floor(0.1*N):size(I,2))));
pSig4 = max(max(I(:,1:floor(0.1*N))));
pSig = max([pSig1 pSig2 pSig3 pSig4]);
numElem = numel(J2(1:floor(0.1*M),:)) + numel(J2(:,1:floor(0.1*N))) + numel(J2(:,size(G,2)-floor(0.1*N):size(G,2))) + numel(J2(size(G,1)-floor(0.1*M):size(G,1),:)); 
noiseSqr = sqrt(sum(sum((J2(1:floor(0.1*M),:)-I(1:floor(0.1*M),:)).^2)));
noiseSqr = noiseSqr + sqrt(sum(sum((J2(:,1:floor(0.1*N))-I(:,1:floor(0.1*N))).^2)));
noiseSqr = noiseSqr + sqrt(sum(sum((J2(:,size(G,2)-floor(0.1*N):size(G,2))-I(:,size(G,2)-floor(0.1*N):size(G,2))).^2)));
noiseSqr = noiseSqr + sqrt(sum(sum((J2(size(G,1)-floor(0.1*M):size(G,1),:)-I(size(G,1)-floor(0.1*M):size(G,1),:)).^2)));
noiseSqr = noiseSqr / numElem;

pSNREdgeTaper = db(double(pSig)/noiseSqr)
hT2_2 = text(T2Pos(1),T2Pos(2) + size(J2,2)+10,strcat('PSNR =', {' '}, num2str(pSNREdgeTaper,'%.2f'), 'dB'),'HorizontalAlignment','center'); %// Place the text


%% Solving for boundary value problem using tile generation method

tic
alpha = 10;
sigma = 2;
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

vTile = bigIm(M/2+1:M/2+2*M,N/2+1:N/2+2*N);
% 
% % 
% % figure(2)
% % subplot(3,1,1)
% % imshow(Aprime,[])
% % subplot(3,1,2)
% % imshow(Bprime,[])
% % subplot(3,1,3)
% % imshow(Cprime,[])
% % 
% % figure(2)
% % imshow(vTile,[])
% 
% J = deconvlucy(G,PSF);
% figure(3)
% imshow(J,[])
%
figure(3)
imshow(vTile,[])
J = deconvlucy(vTile,PSF,numIt);
J3 = J(M/2+1:M/2+M,N/2+1:N/2+N);
toc
figure(1)
subplot(2,3,6)
imshow(J3,[])
title('Restored Image after using tile method')
hT3 = title('Restored Image after using tile method')
T3Pos = round(get(hT3,'Position')); %// Get the position

%% Numerical analysis




numElem = numel(J3(1:floor(0.1*M),:)) + numel(J3(:,1:floor(0.1*N))) + numel(J3(:,size(I,2)-floor(0.1*N):size(I,2))) + numel(J3(size(I,1)-floor(0.1*M):size(I,1),:)); 
noiseSqr = sqrt(sum(sum((J3(1:floor(0.1*M),:)-I(1:floor(0.1*M),:)).^2)));
noiseSqr = noiseSqr + sqrt(sum(sum((J3(:,1:floor(0.1*N))-I(:,1:floor(0.1*N))).^2)));
noiseSqr = noiseSqr + sqrt(sum(sum((J3(:,size(I,2)-floor(0.1*N):size(I,2))-I(:,size(I,2)-floor(0.1*N):size(I,2))).^2)));
noiseSqr = noiseSqr + sqrt(sum(sum((J3(size(I,1)-floor(0.1*M):size(I,1),:)-I(size(I,1)-floor(0.1*M):size(I,1),:)).^2)));
noiseSqr = noiseSqr / numElem;

pSNRTileAlgo = db(double(pSig)/noiseSqr)
hT3_2 = text(T3Pos(1),T3Pos(2) + size(J3,2)+10,strcat('PSNR =', {' '}, num2str(pSNRTileAlgo,'%.2f'), 'dB'),'HorizontalAlignment','center'); %// Place the text


end