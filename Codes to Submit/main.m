function main()

clear all;
clc;

I = imread('thank.jpg');

if numel(size(I)) >= 3
    I = rgb2gray(I);
end

figure(5)
imshow(I)
numIte = 10;
theta1 = 90;    %Horizontal Motion Blur
theta2 = 0;     %Vertical Motion Blur
lengthB = 50;

[I1, I2, dBI1, dBI2] = Assignment_1(I, numIte, lengthB, theta1, theta2);

% figure()
% subplot(2,3,1)
% imshow(I)
% title('Original Image')  
% subplot(2,3,2)
figure(1)
imshow(I1)
% title('Image1 with Horizontal Motion Blur')  
% subplot(2,3,3)
figure(2)
imshow(I2)
% title('Image2 with Vertical Motion Blur')  
% subplot(2,3,5)
figure(3)
imshow(dBI1,[])
% title(strcat('Restored Image1 after', {' '}, num2str(numIte), ' iterations'))  
% subplot(2,3,6)
figure(4)
imshow(dBI2,[])
% title(strcat('Restored Image2 after', {' '}, num2str(numIte), ' iterations'))  

% theta1 = -90:1:90;
% theta2 = 0;
% for i = 1:numel(theta1)
%     [I1, I2, dBI1, dBI2] = Assignment_1(I, numIte, 40, theta1(i), theta2);
%     s1(i) = estimate_sharpness(dBI1);
%     s2(i) = estimate_sharpness(dBI2);
%     s3(i) = max([s1(i) s2(i)]);
% end

% figure()
% % subplot(3,1,1)
% plot(theta1, s1)
% % subplot(3,1,2)
% % plot(theta1, s2)
% % subplot(3,1,3)
% % plot(theta1, s3)



end