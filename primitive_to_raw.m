function U = primitive_to_raw(rho,u,theta,nMax)
%PRIMITIVE_TO_RAW Raw moments U_k of a 1-D Maxwellian up to order nMax.
%
% Central normalized Gaussian moments:
% W_0=1, W_1=0, W_2=1, W_4=3, W_6=15, ...

U = zeros(nMax+1,1);
for k = 0:nMax
    val = 0;
    for r = 0:k
        if mod(r,2)==0
            Wr = gaussian_central_moment(r);
            val = val + nchoosek(k,r)*u^(k-r)*theta^(r/2)*Wr;
        end
    end
    U(k+1) = rho*val;
end
end
