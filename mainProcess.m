clc
clear all

%%  MA process 
N = 3500; %number of data
np = 0.01; %noise power
sp = 1; %signal power
h = [1 2.8 -3 ]; % 3 tap filter - true impulse response
u = sqrt(sp/2).*randn(1,N+1);
x = u(1:N) + u(2:N+1);  % Moving average process :: INPUT
d = conv(x,h);  % True target
d = d(1:N) + sqrt(np).*randn(1,N); % noisy target :: OUTPUT
mu = 0.002;

%% Objective is to get the impulse response :: SYSTEM IDENTIFICATION :: LMS
option = 'Simple';
sz = length(h); %no of taps
mu = 0.003;
[~,Ws,es] = computeLMS(x,d,mu,sz,N,option);


%%
option = 'Ideal';
sz = length(h); %no of taps
mu = 0.002;
[~,Wideal,eideal] = computeLMS(x,d,mu,sz,N,option);

%%
option = 'DCT';
sz = 3; %no of taps
mu = 0.002;
[~,Wdct,edct] = computeLMS(x,d,mu,sz,N,option);

%%
sz = length(h); %no of taps
mu = 0.003;
iter = 100;
[~,Wpcg,epcg, U_opt, cnd] = computePrecoGLMS(x',d,mu,sz,N,iter);
%% Objective is to get the impulse response :: SYSTEM IDENTIFICATION :: LMS




%% _________________________________________________________ %%
%% Autoregressive process :: 1st order MARKOV%%
clc
clear all

N = 3200;
np=0.01;
sp = 1;

% generate X
x = zeros(N,1);
x(1) = sqrt(5)*randn(1);
rho = 0.9;
for pp = 2:N
   x(pp) = sqrt(5)*randn(1)-rho*x(pp-1);
end
x = sqrt(sp)*x;
h = [1 -0.8 6 3];   % These are to be estimated
d = conv(x,h);
d = d(1:N) + sqrt(np).*randn(1,N); % noisy target :: OUTPUT

%%
option = 'Simple';
sz = length(h); %no of taps
mu = 0.004;
[y1,W,e] = computeLMS(x,d,mu,sz,N,option);

%%
option = 'DCT';
sz = 3; %no of taps
mu = 0.004;
[y1dct,W,edct] = computeLMS(x,d,mu,sz,N,option);

%% 
option = 'DCT';
sz = 3; %no of taps
mu = 0.004;
[y1dct,W,edct] = computeLMS(x,d,mu,sz,N,option);


%%
option = 'Ideal';
sz = 3; %no of taps
mu = 0.004;
[y1ideal,w0ideal,w1ideal,w2ideal,eideal] = computeLMS(x,d,mu,sz,N,option);


%% 
sz = length(h); %no of taps
mu = 0.004;
iter = 100;
[y1pcg, Wpcg ,epcg, U_opt, cnd] = computePrecoGLMS(x,d,mu,sz,N,iter);


%%_______________________________________%%
%% 2nd order autoregressive
clc
clear all

%%
N = 2000;
np=0.01;
sp = 1;
rho1 = 0.9;
rho2 = 0.2;

% generate X
x = zeros(N,1);
ss = 5;
x(1) = sqrt(ss)*randn(1);
x(2) = sqrt(ss)*randn(1)-rho1*x(1);

for pp = 3:N
   x(pp) = sqrt(ss)*randn(1)-rho1*x(pp-1)-rho2*x(pp-2);
end
x = sqrt(sp)*x;
h = [1 6 -3];   % This are to be estimated
d = conv(x,h);
d = d(1:N) + sqrt(np).*randn(1,N); % noisy target :: OUTPUT


%%
option = 'Simple';
sz = length(h); %no of taps
mu = 0.0005;
[~,Ws,es] = computeLMS(x,d,mu,sz,N,option);


%%
option = 'DCT';
sz = 3; %no of taps
mu = 0.002;
[~,Wdct,edct] = computeLMS(x,d,mu,sz,N,option);

%%
sz = length(h); %no of taps
mu = 0.002;
iter = 120;
[~,Wpcg,epcg, U_opt, cnd] = computePrecoGLMS(x,d,mu,sz,N,iter);












%% _________________________________%%
%% Hebb-LMS
%% process
clc
clear all
%%
N = 1500;
np=0.01;
sp = 1;

%%
clc
clear all
N = 3200;
np=0.01;
sp = 1;
x1 = sqrt(4)*randn(N/2,1)+1 + sqrt(np).*randn(N/2,1);
x2 = sqrt(4)*randn(N/2,1)-2 - sqrt(np).*randn(N/2,1);

x = [x1;x2];


%% weights
Wi = randn(1,4);
option = 'Simple';
sz=4;
mu=0.001;
W_t = Wi;
e_t_s = [];
for k = 1:500
    xx = x(randperm(length(x)));
    %xx=x;
    [ys,W,e] = computeHebbLMS1D(xx,mu,sz,N,Wi,option);    
    W_t = [W_t;W];
    e_t_s = [e_t_s;e];
    Wi = W;
end

clear W k Wi e
%% PrecoG
Wi = randn(1,4);
status = 'N';
sz=4;
mu=0.001;

U = eye(sz,sz); %for first iter

W_t = Wi;
e_tt = [];
iter = 200;
for kk = 1:500
    xx = x(randperm(length(x)));
    %xx = x;
    [ytt,W,e, U_opt] = computePrecoGHebbLMS1D(xx,mu,sz,N,Wi,U,iter,status);
    W_t = [W_t;W];
    e_tt = [e_tt;e];
    status = 'Y';
    U = U_opt;
    Wi = W;
end
%% _________________________%%

