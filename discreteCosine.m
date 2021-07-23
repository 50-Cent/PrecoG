function [op, U] = discreteCosine(ip)

[~,N] = size(ip);

%% Discrete Cosine transformation matrix
U_dc = zeros(N,N);
for i = 1:N
    for l = 1:N
        if i==1 
            U_dc(i,l) = sqrt(1/N)*cos(((i-1)*(l-1+0.5)*pi)/N);
        else
            U_dc(i,l) = sqrt(2/N)*cos(((i-1)*(l-1+0.5)*pi)/N);
        end
    end
end

op = U_dc*ip*U_dc';
U = U_dc';