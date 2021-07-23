function [op] = discreteFourier(ip)

[~,N] = size(ip);

%% Discrete Cosine transformation matrix
U_df = zeros(N,N);
for k = 1:N
    for l = 1:N         
        U_df(k,l) = sqrt(1/N)*exp(1i*((2*pi*k*l)/N)) ;       
    end
end

op = U_df*ip*U_df';