h = zeros(32,1,10);
h_cap = zeros(32,1,10);
h_cap_c = zeros(32,1,10);
h_cap_c_2 = zeros(32,1,10);
h_cap_5 = zeros(6,1,10);

for t = 1:10
bits = randi([0,1],1024,1);
symbol = zeros(512,1);
for i = 1:512
   symbol(i) = ((1-2*bits(2*i-1))+1j*(1-2*bits(2*i)));
end

X = diag(symbol);
F = zeros(512,32);
for i = 1:512
    for j = 1:32
        F(i,j)=exp(2*pi*(i-1)*(j-1)*1j/512);  
    end
end

h(:,:,t) = zeros(32,1);
p = zeros(32,1);
lambda = 0.2;
std = (0.5)^(0.5);

a = std*randn(32,1);
b = std*randn(32,1);
for i = 1:32 
    p(i) = exp(-1*lambda*(i-1));
end

p_norm=norm(p);
for i = 1:32
    h(i,1,t) = (a(i)+1j*b(i))*p(i)/p_norm;
end
% zero = randperm (32,26);
% for i = 1:26
% h(zero(i),1,t)=0;
% end

noise = randn(512,1)+1j*randn(512,1);

for i = 1:512
noise(i) = noise(i)/abs(noise(i));
end
Y = mtimes(X,F);
y = Y*h(:,:,t) + noise;
P = mtimes(ctranspose(Y),Y);
Q1 = inv(P);
Q = mtimes(Q1,ctranspose(Y));
h_cap(:,:,t) = Q*y;
%plot(h,h_cap);

% when we know the sparsity information %
% A = zeros(32,32);
% for i =1:26
%     A(zero(i),zero(i))=1;
% end
% b_bar = zeros(32,1);
% Z = mtimes(A,Q1);
% l = 2*inv(mtimes(Z,transpose(A)))*(A*h_cap(:,:,t)-b_bar);
% h_cap_c(:,:,t) = h_cap(:,:,t) - 0.5*mtimes(Q1,transpose(A))*l;
% plot(h(27:end),h_cap(27:end));

% question 5 %

%     A = mtimes(X,F);
%     r = y;
%     support = zeros(32,1);
% for k = 1:6
%     [val1, index] = max(mtimes(ctranspose(A),r));
%     support(k) = index;
%     support = unique(support);
%     [val2, i1] = min(support);
%     [val3, i2] = max(support);
%     B = A(:,i1:i2);
%     C = mtimes(inv(mtimes(ctranspose(B),B)),ctranspose(B));
%     P_k = mtimes(B,C);
%     r = (eye(512)-P_k)*y;
%     
% end
% P = mtimes(ctranspose(B),B);
% Q11 = inv(P);
% Q = mtimes(Q11,ctranspose(B));
% %h_cap_5(:,:,t) = Q*y;




% question 4 %
A_c = zeros(3,32);
A_c(1,1)=1; A_c(1,2)=-1;
A_c(2,3)=1; A_c(2,4)=-1;
A_c(3,5)=1; A_c(3,6)=-1;
Z = mtimes(A_c,Q1);
b_bar=zeros(3,1);
l = 2*inv(mtimes(Z,transpose(A_c)))*(A_c*h_cap(:,:,t)-b_bar);
h_cap_c_2(:,:,t) = h_cap(:,:,t) - 0.5*mtimes(Q1,transpose(A_c))*l;
%plot(h,h_cap_c_2);
end
plot(h(:,:,1),h_cap(:,:,1));
E_h=sum(h,3);
E_h_cap=sum(h_cap,3);
E_h_cap_c=sum(h_cap_c,3);
E_h_cap_c_2=sum(h_cap_c_2,3);
plot(E_h,E_h_cap);