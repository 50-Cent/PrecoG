function [y1,W,e] = computeLMS(x,d,mu,sz,N,option)

%%
if isequal(option, 'Ideal')
%     R = [1 0.5 0;0.5 1 0.5;0 0.5 1];  %because we know it is an MA(2) process.
    if sz==3
%         R = [1 0.9 0.81;0.9 1 0.9; 0.81 0.9 1];  % AR(1)
        R = [1 0.5 0;0.5 1 0.5;0 0.5 1]; 
    elseif sz==4
       R = [1 0.9 0.81 0.729;0.9 1 0.9 0.81; 0.81 0.9 1 0.9; 0.729 0.81 0.9 1]; 
    end
    [U,V] = eig(R);
    V = diag(V);    
    U = U';
elseif isequal(option, 'Simple')
    U = eye(sz,sz);
elseif isequal(option, 'DCT')
    U = zeros(sz,sz);
    for i = 1:sz
        for l = 1:sz
            if i==1 
                U(i,l) = sqrt(1/sz)*cos(((i-1)*(l-1+0.5)*pi)/sz);
            else
                U(i,l) = sqrt(2/sz)*cos(((i-1)*(l-1+0.5)*pi)/sz);
            end
        end
    end
elseif isequal(option, 'DFT')
    U = zeros(sz,sz);
    for k = 1:sz
        for l = 1:sz         
            U(k,l) = sqrt(1/sz)*exp(1i*((2*pi*k*l)/sz)) ;       
        end
    end    
end
%%
w0 = zeros(N,1); w1 = zeros(N,1); w2 = zeros(N,1); w3 = zeros(N,1);
y1 = zeros(N,1); e = zeros(N,1); P = zeros(N,1); 
bta = 0.9;
epss = 0.001;


if sz == 3
    P(1) = (1-bta)*(x(1)^2);

    y1(1) = w0(1)*x(1);
    e(1) = d(1)-y1(1);
    w0(2) = w0(1)+2*mu*e(1)*x(1);

    P(2) = bta*P(1) + (1-bta)*(x(2)^2);
    %
    y1(2) = w0(2)*x(2)+w1(2)*x(1);
    e(2) = d(2)-y1(2);
    w0(3) = w0(2)+2*mu*e(2)*x(2);
    w1(3) = w1(2)+2*mu*e(2)*x(1);
elseif sz==4
        P(1) = (1-bta)*(x(1)^2);

        y1(1) = w0(1)*x(1);
        e(1) = d(1)-y1(1);
        w0(2) = w0(1)+2*mu*e(1)*x(1);

        P(2) = bta*P(1) + (1-bta)*(x(2)^2);
        %
        y1(2) = w0(2)*x(2)+w1(2)*x(1);
        e(2) = d(2)-y1(2);
        w0(3) = w0(2)+2*mu*e(2)*x(2);
        w1(3) = w1(2)+2*mu*e(2)*x(1);
        P(3) = bta*P(2)+(1-bta)*(x(3)^2);
        
        y1(3) = w0(3)*x(3)+w1(3)*x(2)+w2(3)*x(1);
        e(3) = d(3)-y1(3);
        w0(4) = w0(3)+2*mu*e(3)*x(3);
        w1(4) = w1(3)+2*mu*e(3)*x(2);
        w2(4) = w2(3)+2*mu*e(3)*x(1);      
end
    
if sz==3
    for ii = 3:N-1

            tmp = U*[x(ii) x(ii-1) x(ii-2)]'; % transform
            if ~isequal(option, 'Simple')
                P(ii) = bta*P(ii-1) + (1-bta)*(tmp(1)^2); % power of transformed x(ii)    
                vv = tmp./[sqrt(P(ii)+epss) sqrt(P(ii-1)+epss) sqrt(P(ii-2)+epss)]'; %normalization  
            else 
                vv = tmp;
            end
            wgt = [w0(ii) w1(ii) w2(ii)];       

            y1(ii) = wgt*vv;    
            e(ii) = d(ii)-y1(ii);    

            w0(ii+1) = w0(ii)+2*mu*e(ii)*vv(1);
            w1(ii+1) = w1(ii)+2*mu*e(ii)*vv(2);
            w2(ii+1) = w2(ii)+2*mu*e(ii)*vv(3);   
    end
    
elseif sz==4
    for ii = 4:N-1
            tmp = U*[x(ii) x(ii-1) x(ii-2) x(ii-3)]'; % transform
            P(ii) = bta*P(ii-1) + (1-bta)*(tmp(1)^2);
            vv = tmp./[sqrt(P(ii)+epss) sqrt(P(ii-1)+epss) sqrt(P(ii-2)+epss) sqrt(P(ii-3)+epss)]'; %normalization  
            wgt = [w0(ii) w1(ii) w2(ii) w3(ii)];    
            y1(ii) = wgt*vv;    
            e(ii) = d(ii)-y1(ii);    

            w0(ii+1) = w0(ii)+2*mu*e(ii)*vv(1);
            w1(ii+1) = w1(ii)+2*mu*e(ii)*vv(2);
            w2(ii+1) = w2(ii)+2*mu*e(ii)*vv(3); 
            w3(ii+1) = w3(ii)+2*mu*e(ii)*vv(4); 
     end
    
    
end
if sz == 3
    W = [w0 w1 w2];
elseif sz == 4
    W = [w0 w1 w2 w3];
end

e = abs(e);
end