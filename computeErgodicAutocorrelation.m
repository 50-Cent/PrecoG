function M = computeErgodicAutocorrelation(x,sz)

M = zeros(sz,sz);
R = zeros(sz,1);

mu = mean(x);
sd = std(x);

x = (x-mu)/sd;
length(x)

for kk = 1:sz
   tmp = [x(kk:length(x)); zeros(kk-1,1)]; 
   length(tmp)
   lnn = length(x)-kk+1;
   R(kk) = ((x')*tmp)/lnn;
end


for ii = 1:sz
    for jj = ii:sz
        M(ii,jj) = R(jj-ii+1);
        M(jj,ii) = M(ii,jj);
    end
end


end