% Question 1 %

mu = zeros(2,1);
sigma = [[1,0];[0,2]];
n = 100;
R = mvnrnd(mu,sigma,n);
R = transpose(R);
% MLE
mean = sum(R,2)/n;

X = zeros(2,2);
for i=1:n
    
   X = X + (R(:,i)-mean)*transpose(R(:,i)-mean);
    
end
covariance = X/n

% Question 4 %

% Inv Wishart Distribution %

df = 5;
Tau = [[2,0];[0,4]];
m=100000;
prior = zeros(2,2,m);
for i=1:m
    
    prior(:,:,i) = iwishrnd(Tau,df);

end

numerator = zeros(2,2);
denominator = 0;
P = zeros(m,1);

for i=1:m
    
    P(i) = sum(sum(transpose(R)*inv(prior(:,:,i)).*transpose(R)));
    denominator = denominator + det(prior(:,:,i))^(-n/2)*exp(-0.5*P(i))/m
    numerator = numerator + prior(:,:,i)*det(prior(:,:,i))^(-n/2)*exp(-0.5*P(i))/m
    
end

A = numerator/denominator
