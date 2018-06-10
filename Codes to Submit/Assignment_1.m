function [I1, I2, dBI1, dBI2] = Assignment_1(I, numIte, len, theta1, theta2)

%Manually Introducing Horizontal Motion Blur to Image I
PSFV = fspecial('motion', len, theta1);
I1 = imfilter(I, PSFV, 'conv', 'circular');

%Manually Introducing Vertical Motion Blur to Image I
PSFH = fspecial('motion', len, theta2);
I2 = imfilter(I, PSFH, 'conv', 'circular');

G1 = fft2(I1);
G1C = conj(G1);
G2 = fft2(I2);
G2C = conj(G2);

%Iterating over equations 11 to 14 for restoration
for ite = 1:numIte    
    
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
    
    G1 = G11;
    G2 = G21;
    G1C = conj(G1);
    G2C = conj(G2);
end

dBI1 = ifft2(G1); %deblurred Image from I1
dBI2 = ifft2(G2); %deblurred Image from I2
end