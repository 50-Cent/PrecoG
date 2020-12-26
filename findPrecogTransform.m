function [min_cnd, U_opt] = findPrecogTransform(R,A,epsn1,epsn2,l2coeff,learning_rate,maxIter,log_penalty)

[indX, indY] = find(triu(A)~=0); %only the upper triangle (A is symmetric)
%maxIter = 200;
N = size(R,1);
stp = 0.001;
%l2coeff = 0.18;

cnd_num = [];
[~,bb] = eig(R);
bb = diag(bb);
cnd_num = [cnd_num;max(bb)/min(bb)];
min_cnd = 99999;

clear aa bb

for i = 1:maxIter
   D = diag(sum(A,1));
   L = D-A;
   L = inv(sqrt(D))*L*inv(sqrt(D));   %single connected
   [U,Dg] = eig(L);
   Dg = diag(Dg);
   [Dg,od] = sort(Dg,'ascend');
   U = U(:,od);  %eigenvectors sorted order
     
   clear od
   %% condition number computation
    R_t = U'*R*U;
    pwr = inv(sqrt(diag(diag(R_t))));
    R_t = pwr*R_t*pwr;
    [~,Eg] = eig(R_t);
    Eg = diag(Eg);
    cnd_num = [cnd_num;max(Eg)/min(Eg)];
    if isreal(cnd_num(end))
        if cnd_num(end)>= 1
            if cnd_num(end)<min_cnd
                U_opt = U;
                min_cnd = cnd_num(end);
            end
        end
    end
   
   %% truncation of eigenvalues
   Dg = (floor(Dg*1000))/1000; %MATLAB problem
   
   %% calculation  of dE/dU   
   I = eye(size(R,1));
   S = U'*R;
   dE_dU = 2*( (2*R - 2*(2+epsn1-epsn2)*I + ( (1+epsn1)^2+(1+epsn2)^2) )*((S*U).*I) )*S';
      

   %% Clustering of eigenvalues :: 
   % Dg is sorted. We can take advantage of that.
   % 0: for zero eigenvalue, 1: eigenvalue with multiplicity 1, other: eigenvalue w multiplicity >1
   
   eigStatus = eigenvalueGroup(Dg);
   loc0 = find(Dg==0);
   eigStatus(loc0) = 0;  %zero eigenvalues
  
   
   %% du/dw  :: <u1, BdW/dwB^T*u2> for all weights and eigenvectors

   dU_dW = zeros(N,N,length(indX)); %N*N inner product for each w
   for m = 1:length(indX)
      xx = indX(m);
      yy = indY(m);      
      for s = 1:N %eigenvectors
          stt = eigStatus(s);
          if stt ==1  %multiplicity 1           
              lmbd = Dg(s)-Dg;
              lmbd(s)=[];
              lmbd = diag(1./lmbd);
              
              inprod = U(xx,s)*(U(xx,:)-U(yy,:)) + U(yy,s)*(U(yy,:)-U(xx,:));
              inprod(s) = [];
              inprod = diag(inprod);
              
              U_red = U;
              U_red(:,s)=[];
              
              dUs_dWm = sum((U_red*inprod)*lmbd, 2);
              dU_dW(:,s,m) = dUs_dWm;
          elseif stt==0
              loc0 = find(eigStatus==0);
              if length(loc0)==1
                  lmbd = Dg(s)-Dg;
                  lmbd(s)=[];                  
                  lmbd = diag(1./lmbd);
                  
                  inprod = U(xx,s)*(U(xx,:)-U(yy,:)) + U(yy,s)*(U(yy,:)-U(xx,:));
                  inprod(s) = [];
                  inprod = diag(inprod);
                  
                  U_red = U;
                  U_red(:,s)=[];
                  
                  dUs_dWm = sum((U_red*inprod)*lmbd, 2);
                  dU_dW(:,s,m) = dUs_dWm;                 
              else
                  remLoc0 = setdiff(loc0,s); %if more than one 0 eigenvalue 
                  lmbd = Dg(s)-Dg;                   
                  lmbd(remLoc0) = 1;
                  lmbd(s)=[];
                  lmbd = diag(1./lmbd); %N-1 x N-1
                  
                  inprod = U(xx,s)*(U(xx,:)-U(yy,:)) + U(yy,s)*(U(yy,:)-U(xx,:));
                  inprod(s) = [];
                  inprod = diag(inprod); %N-1 x N-1
                  
                  U_red = U;                  
                  eig2replace = repmat(U_red(:,s),1,length(remLoc0));
                  U_red(:,remLoc0) = eig2replace;
                  U_red(:,s)=[]; % N x N-1
                  
                  dUs_dWm = sum((U_red*inprod)*lmbd, 2);
                  dU_dW(:,s,m) = dUs_dWm;                              
              end
          else
              loc0 = find(eigStatus==stt);
              remLoc0 = setdiff(loc0,s); %if more than one 0 eigenvalue 
              lmbd = Dg(s)-Dg;              
              lmbd(remLoc0) = 1/Dg(s);
              lmbd(s)=[];
              lmbd = diag(1./lmbd); %N-1 x N-1

              inprod = U(xx,s)*(U(xx,:)-U(yy,:)) + U(yy,s)*(U(yy,:)-U(xx,:));
              inprod(s) = [];
              inprod = diag(inprod); %N-1 x N-1

              U_red = U;                  
              eig2replace = repmat(U_red(:,s),1,length(remLoc0));
              U_red(:,remLoc0) = eig2replace;
              U_red(:,s)=[]; % N x N-1

              dUs_dWm = sum((U_red*inprod)*lmbd, 2);
              dU_dW(:,s,m) = dUs_dWm;         
              
          end
          
          % write update           
          hh = A(xx,yy);
          A(xx,yy) = A(xx,yy)*(1-2*l2coeff)- learning_rate*trace( (dE_dU)'*dU_dW(:,:,m) );
          if A(xx,yy)<0
              A(xx,yy) = A(xx,yy)+ log_penalty*(1/hh);  %log compensation
          end
              
          A(yy,xx) = A(xx,yy);
          
      end  
   end  
    
end
% cnd_num = real(cnd_num);
% cnd_num = cnd_num(cnd_num>1);
% cnd_num = min(cnd_num);
%min(cnd_num)
%disp('... Done...\n')
end

%% Examples found experimentally
% size 3 :: MARKOV-1
% sweet spot =  log compen=0.01, learning=0.02, L2 compen = 0.1, A =
% structured

%size 3 :: MARKOV -I (rho=0.15)
% cnd_DC=1.12, cnd_PrecoG = 1.02 (learning = 0.03, l2 compen = 0.42, iter = 10)


%size 3 :: MARKOV -I (rho=0.05)
% cnd_DC=1.15, cnd_PrecoG = 1.04 (learning = 0.0035, l2 compen = 0.2, iter=100)

%size 3 :: MARKOV-II (0.9,0.7)
% cnd_DC = 1.23, cnd_PrecoG = [1.08,1.05] (learning = [0.0025, 0.00034], l2 compen = 0.2, iter=100)

% size 5: HILBERT 
%  cnd_DC=41.46, cnd_PrecoG = [4.17,7.43] (learning = [0.001,0.003], l2 compen = 0.18, iter=200)

% size 4: MARKOV-I (0.96)
% cnd_DC= 1.0012, cnd_PrecoG = 1.0000 (grid search)

% % size 4: HILBERT (0.001) 
%  cnd_DC=107.67, cnd_PrecoG = [8.24] (grid search)

