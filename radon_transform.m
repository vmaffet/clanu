function I_radon = radon_transform(Im,N_t,N_theta)

I_radon = zeros(N_theta,N_t);
[N1,N2] = size(Im);

Theta = linspace(0,pi,N_theta+1);
T = linspace(-1,1,N_t);

Nd = max(N1,N2);
delta_t = linspace(-1,1,Nd);

for i_theta = 1:N_theta,
    
    for j_t = 1:N_t,
        
        ppx_theta =   cos(Theta(i_theta))*delta_t - sin(Theta(i_theta))*T(j_t);
        ppy_theta =   sin(Theta(i_theta))*delta_t + cos(Theta(i_theta))*T(j_t);
        sum=0;
        
        for k=1:Nd
          if (ppx_theta(k)^2 + ppy_theta(k)^2<0.98),
       
       indice_x = (ppx_theta(k)+1)*(N1-1)/2+1;
       indice_y = (ppy_theta(k)+1)*(N2-1)/2+1;
       lambda_x = indice_x -floor(indice_x);
       lambda_y = indice_y -floor(indice_y);
       
       sum=sum + (1-lambda_x)*(1-lambda_y)*Im(floor(indice_x),floor(indice_y)) ...
           +(lambda_x)*(1-lambda_y)*Im(min(floor(indice_x)+1,N1),floor(indice_y)) ...
           +(1-lambda_x)*(lambda_y)*Im(floor(indice_x),min(floor(indice_y)+1,N2)) ...
           +(lambda_x)*(lambda_y)*Im(min(floor(indice_x)+1,N1),min(floor(indice_y)+1,N2)) ; 
           end
        end
        
         I_radon(i_theta,j_t) = sum/Nd;
    
    end
end

end