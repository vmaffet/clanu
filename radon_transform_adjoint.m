function Im = radon_transform_adjoint(I_radon,N1,N2)

Im =zeros(N1,N2);
[X2,X1] = meshgrid(linspace(-1,1,N1),linspace(-1,1,N2));
[N_theta,N_t] = size(I_radon);
Theta = linspace(0,pi,N_theta +1);
T = linspace(-1,1,N_t);

for i_theta = 1:N_theta,

     t_theta =  X2*cos(Theta(i_theta)) - X1*sin(Theta(i_theta));
     
     indice_t = (t_theta+1)/2*(N_t-1)+1;
     lambda_t = indice_t - floor(indice_t);
           
     for i=1:N1,
         for j=1:N2,
             if  (sqrt(X1(i,j).^2 + X2(i,j).^2)<0.98) 
             if (t_theta(i,j)>-0.98 &&  t_theta(i,j)<0.98),
                Im(i,j) =  Im(i,j)  + ((1-lambda_t(i,j))*I_radon(i_theta,floor(indice_t(i,j))) ...
                               +  lambda_t(i,j)*I_radon(i_theta,min(floor(indice_t(i,j))+1,N_t)))*2*pi/N_theta;
             end
             end
         end
     end
        
end