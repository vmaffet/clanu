% Flavien Chenost
% Vincent Maffet
clear all;
close all;

N1= 2^7; % y output pixels
N2= N1;  % x output pixels 
Nt= N1;  % number of sensors
No= 30;  % number of angle steps
Ni= 50;  % number of iterations

% TEST DATA %
% a= ones(N1/8,1);
% b= zeros(N1/8,1);
% x= [b' b' b' a' a' b' b' b'];
% y= [b' b' b' b' a' b' b' b'];
% [X Y]= meshgrid(x, y);
% I= X.*Y;

% REAL APPLICATION %
I= phantom(N1, N2);

% DATA AQUISITION %
Idata= radon_transform(I, Nt, No);
Idata= Idata + 0.01*randn(No, Nt); % adding some noise %

% Optionnal Filtering %
k= 1/2*[0:Nt/2,-Nt/2+1:-1];
% for i= 1:No
%     Idata(i,:)= ifft((16*cos(4*pi*k/Nt)+16).*fft(Idata(i,:)));
% end
% 32-abs(k)
% 32*cos(2*pi*k/Nt)
% 16*cos(4*pi*k/Nt)+16 LE MEILLEUR

% RECONSTRUCTION %
tau= 0.3; % Limite= 0.38 pour N= 128
A= zeros(N1, N2);
B= zeros(N1, N2);
for i= 1:Ni
    A= A - tau*radon_transform_adjoint(radon_transform(A, Nt, No) - Idata, N1, N2);
    B= B - tau*radon_transform_adjoint(radon_transform(B, Nt, No) - Idata, N1, N2);
    for k= 1:N1
        for l= 1:N2
            if (A(k,l) < 0) 
                A(k,l)= 0;
            end
        end
    end
end

% Display %
subplot(2,2,1)
imagesc(I);
title("Input body");
subplot(2,2,2)
imagesc(Idata);
title("Acquired data");
subplot(2,2,3)
imagesc(B);
title("Vanilla Output");
subplot(2,2,4)
imagesc(A);
title("Forced positive");