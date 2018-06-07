function redBoundEff()

I = imread('cameraman.tif');

if numel(size(I)) >= 3
    I = rgb2gray(I);
end

%Create a sample image with noise.

% Set the random number generator back to its default settings for
% consistency in results.
rng default;

% I = checkerboard(8);
PSF = fspecial('gaussian',10,10);
V = .0000;
BlurredNoisy = imnoise(imfilter((I),PSF),'gaussian',0,V);

% Create a weight array to specify which pixels are included in processing.

WT = zeros(size(I));
WT(5:end-4,5:end-4) = 1;
INITPSF = ones(size(PSF));

% Perform blind deconvolution.

[J P] = deconvblind(double(BlurredNoisy),INITPSF,20,10*sqrt(V),WT);

figure()
imshow(I);

figure()
imshow(J,[]);
end