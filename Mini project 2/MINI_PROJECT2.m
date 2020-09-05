%%
N = [1,10,100,1000,10000];
A = 10;
E_A = zeros([5,1]);
iterations = 5000;
for j= 1:5
    U = zeros([N(j),iterations]);
    n = zeros([N(j),iterations]);
    x = zeros([N(j),iterations]);
    for i= 1:iterations
        %generating Laplace distribution
        U(:,i) = rand([N(j),1]) -0.5;
        n(:,i) = -1.*sign(U(:,i)).*log(1-2.*abs(U(:,i)))./(sqrt(2));
        x(:,i) = A + n(:,i);
        A_cap(i,j) = median(x(:,i));
        E_A(j) = E_A(j) + A_cap(i,j);
    end
    E_A(j) = E_A(j)/iterations;
    var(j) = (sum((A_cap(:,j)-E_A(j).*ones([iterations,1])).^2))/iterations;
end
%%
%% PLOTS %%
%%
x_axis = [1:1:5];
stem(x_axis,E_A,'r','LineWidth',2);
hold on
stem(x_axis,A*ones([5,1]),'b--','LineWidth',1);
ylim(A*[0.9,1.1]);
xlabel('Number of samples N = [1,10,100,1000,10000]');
ylabel('Expectation of MLE of A and A');
title('Comparison of expectation of MLE of A and A');
legend('Expectation of MLE of A', 'A');
%%

stem(x_axis,var,'r','LineWidth',2);
xlabel('Number of samples N = [1,10,100,1000,10000]');
ylabel('Variance of MLE of A');
title('Variance vs number of samples');

%%
[c,x_axis_c] = ecdf(A_cap(:,5));
[c2,x_axis_c] = ecdf(sqrt(10000)*(A_cap(:,5)-E_A(5)*ones([iterations,1])));
gauss_fit = normpdf(x_axis_c,0,1/2);
[c3,x_axis_c] = ecdf(gauss_fit);
plot(x_axis_c,c,'g');
xlabel('A');
ylabel('CDF(A-cap)');
title('CDF of the estimate of A');
%%

plot(x_axis_c,c2,'b');
hold on
plot(x_axis_c,c3, 'r--');
xlabel('x');
ylabel('CDF(x)');
title('Normal convergence of MLE in distribution');
legend('CDF(root(N)*(A-cap-E[A]))', 'CDF(N(0,1/I(A)))');
%%

[pdf_A,x_axis_c]=hist(A_cap(:,5));
plot(x_axis_c,pdf_A/N(5),'r');
xlabel('A cap');
ylabel('PDF(A cap)');
title('PDF of the MLE of A');
%%

[pdf,x_axis_c]=hist(sqrt(N(5))*(A_cap(:,5)-E_A(5)*ones([iterations,1])));

plot(x_axis_c,pdf/iterations,'r');
hold on
gauss_fit = normpdf(x_axis_c,0,1/2);
plot(x_axis_c,gauss_fit,'b');
xlabel('x');
ylabel('PDF(x)');
title('Normal convergence of MLE in distribution');
legend('PDF(root(N)*(A-cap-E[A]))', 'PDF(N(0,1/I(A)))');