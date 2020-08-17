clc
clear
X1=csvread('.csv');
X2=csvread('.csv');
X=[X1,X2];
M=mapminmax(X,0,1);
V=M;
[m n]=size(V);  %计算V的规格  
r=270;
W=rand(m,r);  %初始化W H
W=W./(ones(m,1)*sum(W));  %归一化W的每一列
H=rand(r,n);
maxiter=5000;  %最大迭代次数 
ONES=ones(m,n);
for iter=1:maxiter 
    W = W.*((V./(W*H))*H')./(ONES*H');
    W=W./(ones(m,1)*sum(W)); 
    H = H.*(W'*( V./(W*H)))./(W'*ONES); 
end
csvwrite('\.csv',W);       
  
 

