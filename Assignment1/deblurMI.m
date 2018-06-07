function deblurMI()

% I = imread('DSC_0619.jpg');
I = imread('peppers.png');
if numel(size(I)) >= 3
    I = rgb2gray(I);
end
sharpness=estimate_sharpness(double(I));
% disp(['Sharpness of original image: ' num2str(sharpness)]);
ite = 1;

% for THETA = -90:1:90
LEN = 10;
THETA = 90;
PSFH = fspecial('motion', LEN, THETA);
I2 = imfilter(I, PSFH, 'conv', 'circular');
sharpness=estimate_sharpness(double(I2));
% disp(['Sharpness of Vertical Blur: ' num2str(sharpness)]);
figure()
imshow(I2);
% title('Ver Blurred Image');

% for THETA = -90:1:90

THETA2 = 40;


%     ite
%     
LEN2 = 10;
% THETA2 = 90;
PSFV = fspecial('motion', LEN2, THETA2);
I1 = imfilter(I, PSFV, 'conv', 'circular');
% sharpness=estimate_sharpness(double(I1));
% disp(['Sharpness of Horizontal Blur: ' num2str(sharpness)]);
figure()
imshow(I1)
J1 = I2;
J2 = I1;


% theta = 45;
% J1 = imrotate(I2,theta);
% % size(J2)
% % size(I2)
% % J2 = imresize(J2,size(I2));
% J2 = zeros(size(J1));
% size(I2)
% numel(126:509)
% numel(62:573)
% J2(126:509,62:573) = I1;

% G1 = zeros(size(J1));
% G1(126:509,62:573) = fft2(J2);
G1 = fft2(J2);
sg1 = size(G1);
G1C = conj(G1);
G2 = fft2(J1);
% sg2 = size(G2)
G2C = conj(G2);

for i = 1:size(G1,1)
    mult1(i) = 0;
    mult2(i) = 0;
    for k = 1:size(G1,2)
        mult1(i) = mult1(i) + (G2(i,k)*G1C(i,k));
        mult2(i) = mult2(i) + (abs(G1(i,k)))^2;
    end
    
    for j = 1:size(G1,2)
        G11(i,j) = G1(i,j)*mult1(i)./mult2(i);
    end
end


for j = 1:size(G2,2)
    
    mult1(j) = 0;
    mult2(j) = 0;
    for k = 1:size(G2,1)
        mult1(j) = mult1(j) + (G11(k,j)*G2C(k,j));
        mult2(j) = mult2(j) + ((abs(G2(k,j)))^2);
    end
    
    for i = 1:size(G2,1)
        G21(i,j) = G2(i,j)*mult1(j)/mult2(j);
    end
end

G11C = conj(G11);
G21C = conj(G21);

% for i = 1:size(G11,1)
%     mult1(i) = 0;
%     mult2(i) = 0;
%     for k = 1:size(G11,2)
%         mult1(i) = mult1(i) + (G21(i,k)*G11C(i,k));
%         mult2(i) = mult2(i) + (abs(G11(i,k)))^2;
%     end
%     
%     for j = 1:size(G11,2)
%         G12(i,j) = G11(i,j)*mult1(i)./mult2(i);
%     end
% end
% 
% 
% for j = 1:size(G21,2)
%     
%     mult1(j) = 0;
%     mult2(j) = 0;
%     for k = 1:size(G21,1)
%         mult1(j) = mult1(j) + (G12(k,j)*G21C(k,j));
%         mult2(j) = mult2(j) + ((abs(G21(k,j)))^2);
%     end
%     
%     for i = 1:size(G21,1)
%         G22(i,j) = G21(i,j)*mult1(j)/mult2(j);
%     end
% end
% 
% G12C = conj(G12);
% G22C = conj(G22);
% 
% for i = 1:size(G12,1)
%     mult1(i) = 0;
%     mult2(i) = 0;
%     for k = 1:size(G12,2)
%         mult1(i) = mult1(i) + (G22(i,k)*G12C(i,k));
%         mult2(i) = mult2(i) + (abs(G12(i,k)))^2;
%     end
%     
%     for j = 1:size(G12,2)
%         G13(i,j) = G12(i,j)*mult1(i)./mult2(i);
%     end
% end
% 
% 
% for j = 1:size(G22,2)
%     
%     mult1(j) = 0;
%     mult2(j) = 0;
%     for k = 1:size(G22,1)
%         mult1(j) = mult1(j) + (G13(k,j)*G22C(k,j));
%         mult2(j) = mult2(j) + ((abs(G22(k,j)))^2);
%     end
%     
%     for i = 1:size(G22,1)
%         G23(i,j) = G22(i,j)*mult1(j)/mult2(j);
%     end
% end
% 
% G13C = conj(G13);
% G23C = conj(G23);
% 
% for i = 1:size(G13,1)
%     mult1(i) = 0;
%     mult2(i) = 0;
%     for k = 1:size(G13,2)
%         mult1(i) = mult1(i) + (G23(i,k)*G13C(i,k));
%         mult2(i) = mult2(i) + (abs(G13(i,k)))^2;
%     end
%     
%     for j = 1:size(G13,2)
%         G14(i,j) = G13(i,j)*mult1(i)./mult2(i);
%     end
% end
% 
% 
% for j = 1:size(G23,2)
%     
%     mult1(j) = 0;
%     mult2(j) = 0;
%     for k = 1:size(G23,1)
%         mult1(j) = mult1(j) + (G14(k,j)*G23C(k,j));
%         mult2(j) = mult2(j) + ((abs(G23(k,j)))^2);
%     end
%     
%     for i = 1:size(G23,1)
%         G24(i,j) = G23(i,j)*mult1(j)/mult2(j);
%     end
% end
% 
% G14C = conj(G14);
% G24C = conj(G24);
% 
% for iter=1:5
%     for i = 1:size(G14,1)
%         mult1(i) = 0;
%         mult2(i) = 0;
%         for k = 1:size(G14,2)
%             mult1(i) = mult1(i) + (G24(i,k)*G14C(i,k));
%             mult2(i) = mult2(i) + (abs(G14(i,k)))^2;
%         end
% 
%         for j = 1:size(G14,2)
%             G15(i,j) = G14(i,j)*mult1(i)./mult2(i);
%         end
%     end
% 
%     for j = 1:size(G24,2)
% 
%         mult1(j) = 0;
%         mult2(j) = 0;
%         for k = 1:size(G24,1)
%             mult1(j) = mult1(j) + (G15(k,j)*G24C(k,j));
%             mult2(j) = mult2(j) + ((abs(G24(k,j)))^2);
%         end
% 
%         for i = 1:size(G24,1)
%             G25(i,j) = G24(i,j)*mult1(j)/mult2(j);
%         end
%     end
% 
%     G14 = G15;
%     G24 = G25;
%     G14C = conj(G14);
%     G24C = conj(G24);
% 
% end

dBI1It1 = ifft2(G11); 
% s1(ite)=abs(estimate_sharpness(double(dBI1It1)));
% disp(['Sharpness of Horizontal deBlur: ' num2str(sharpness)]);

dBI2It1 = ifft2(G21);
% s2(ite)=abs(estimate_sharpness(double(dBI2It1)));

% sharpness=abs(estimate_sharpness(double(dBI2It1)));
% disp(['Sharpness of Vertical deBlur: ' num2str(sharpness)]);
ite = ite + 1;
% end

% figure()
% imshow(I1)
figure()
% plot(s1)
% figure()
% plot(s2)
imshow(dBI1It1,[])
% % figure(3)
% % % imshow(ifft2(G12),[])
% % % figure(4)
% % % imshow(ifft2(G13),[])
% % % figure(5)
% % % imshow(abs(ifft2(G14)),[])
% figure()
% imshow(I2)
figure()
imshow(dBI2It1,[])
% % % figure(8)
% % % imshow(abs(ifft2(G22)),[])
% % % figure(9)
% % % imshow(abs(ifft2(G23)),[])
% % % figure(10)
% % % imshow(abs(ifft2(G24)),[])
% % ite = ite + 1;
end

% plot(angI,s1,'o-')
% hold on
% plot(angI,s2,'x-')
% hold on
% end