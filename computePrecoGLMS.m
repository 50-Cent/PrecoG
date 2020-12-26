function [y1,W,e, U_opt, cnd] = computePrecoGLMS(x,d,mu,sz,N,iter)

%%
U = eye(sz,sz); %for first iter
w0 = zeros(N,1); w1 = zeros(N,1); w2 = zeros(N,1); w3 = zeros(N,1);
y1 = zeros(N,1); e = zeros(N,1); P = zeros(N,1); 
bta = 0.9;
epss = 0.001;


if sz == 3
    %%
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


    %% first half
    for ii = 3:iter
        tmp = U*[x(ii) x(ii-1) x(ii-2)]'; % transform
        P(ii) = bta*P(ii-1) + (1-bta)*(tmp(1)^2);
        vv = tmp./[sqrt(P(ii)+epss) sqrt(P(ii-1)+epss) sqrt(P(ii-2)+epss)]'; %normalization  
        wgt = [w0(ii) w1(ii) w2(ii)];    
        y1(ii) = wgt*vv;    
        e(ii) = d(ii)-y1(ii);    
    
        w0(ii+1) = w0(ii)+2*mu*e(ii)*vv(1);
        w1(ii+1) = w1(ii)+2*mu*e(ii)*vv(2);
        w2(ii+1) = w2(ii)+2*mu*e(ii)*vv(3);   
    end
    
elseif sz == 4
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
        
        
        %% first half
        for ii = 4:iter
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
%% find ergodic autocorrelation estimate 
R_est = computeErgodicAutocorrelation(x(1:iter),sz);
R_est = abs(R_est);
%% A customized
A = zeros(sz);
for i = 1:sz-1
    for j = i+1:sz       
       A(i,j) = 0.1*abs(i-j);
       A(j,i) = A(i,j);
    end
end

%% find U
epsn1 = 0.2;
epsn2 = 0.2;
log_penalty = 0.01;
maxIter = 125;
U_opt = [];
cnd = 99999;
for l2coeff = 0:0.01:0.3
   for learning_rate = 0.0010:0.0001:0.01 
        [tmpp,U] = findPrecogTransform(R_est,A,epsn1,epsn2,l2coeff,learning_rate,maxIter,log_penalty);
        if tmpp <cnd
            U_opt = U;
            cnd = tmpp;
        end
   end   
end

cnd
U_opt
%%
%% second half

if sz==3
    for ii = iter+1:N-1
        tmp = U_opt*[x(ii) x(ii-1) x(ii-2)]'; % transform
        P(ii) = bta*P(ii-1) + (1-bta)*(tmp(1)^2);
        vv = tmp./[sqrt(P(ii)+epss) sqrt(P(ii-1)+epss) sqrt(P(ii-2)+epss)]'; %normalization  
        wgt = [w0(ii) w1(ii) w2(ii)];    
        y1(ii) = wgt*vv;    
        e(ii) = d(ii)-y1(ii);    

        w0(ii+1) = w0(ii)+2*mu*e(ii)*vv(1);
        w1(ii+1) = w1(ii)+2*mu*e(ii)*vv(2);
        w2(ii+1) = w2(ii)+2*mu*e(ii)*vv(3);   
    end
elseif sz==4
   for ii = iter+1: N-1
       tmp = U_opt*[x(ii) x(ii-1) x(ii-2) x(ii-3)]'; % transform
       P(ii) = bta*P(ii-1) + (1-bta)*(tmp(1)^2);
        vv = tmp./[sqrt(P(ii)+epss) sqrt(P(ii-1)+epss) sqrt(P(ii-2)+epss) sqrt(P(ii-3)+epss)]'; %normalization  
        wgt = [w0(ii) w1(ii) w2(ii) w3(ii)];    
        y1(ii) = wgt*vv;    
        e(ii) = d(ii)-y1(ii);
        
        w0(ii+1) = w0(ii)+2*mu*e(ii)*vv(1);
        w1(ii+1) = w1(ii)+2*mu*e(ii)*vv(2);
        w2(ii+1) = w2(ii)+2*mu*e(ii)*vv(3); 
        w3(ii+1) = w3(ii) + 2*mu*e(ii)*vv(4);
       
   end    
end
  
if sz == 3
    W = [w0 w1 w2];
elseif sz == 4
    W = [w0 w1 w2 w3];
end

e = abs(e);


end
