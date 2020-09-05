%% Initializations
pi=[0.5,0.25];
p=0.35;
q=0.6;
m=[1,10];
n=[10,1000,10000];
%% Q1 Generating Observations
headcount=[];
ntosses=m(2);
ntrials=n(1);
pic=pi(2);
for i = 1:ntrials
    c=rand(1);
    noh=0;
    if c<=pic
        cp=p;
    else
        cp=q;
    end        
    for j = 1:ntosses
        t=rand(1);
        if(t<=cp)
           noh=noh+1; 
        end
    end
    headcount=[headcount,noh];
end

%% Q1 EM Algorithm
est_q1=zeros(1000,2);
iterations=1;
est_q1(iterations,:)=[0.45,0.5];
pi_assume=0.25;
diff=100;
while(diff>=1e-6)
    % E-Step
    c1=[0,0]; c2=[0,0];
    for i = 1:length(headcount)
        t1=pi_assume*(est_q1(iterations,1)^headcount(i)*(1-est_q1(iterations,1))^(ntosses-headcount(i)));
        t2=(1-pi_assume)*(est_q1(iterations,2)^headcount(i)*(1-est_q1(iterations,2))^(ntosses-headcount(i)));
        frac=t1/(t1+t2);
        c1(1)=c1(1)+frac*headcount(i);
        c1(2)=c1(2)+frac*(ntosses-headcount(i));
        c2(1)=c2(1)+(1-frac)*headcount(i);
        c2(2)=c2(2)+(1-frac)*(ntosses-headcount(i));
    end

    % M-Step
    new1=c1(1)/(c1(1)+c1(2));
    new2=c2(1)/(c2(1)+c2(2));
    
    diff=max(abs(new1-est_q1(iterations,1)),abs(new2-est_q1(iterations,2)));
    iterations=iterations+1;
    est_q1(iterations,:)=[new1,new2];
end
%% PLOTS 
stem([1:1:iterations-1],est_q1(1:iterations-1,1));
xlabel('Number of iterations');
ylabel('Estimate at kth iteration');
title('Learning curve of estimate of p');
%%
stem([1:1:iterations-1],est_q1(1:iterations-1,2));
xlabel('Number of iterations');
ylabel('Estimate at kth iteration');
title('Learning curve of estimate of q');


%% Q2 Generating Observations
headcount=[];
ntosses=m(1);
ntrials=n(3);
pic=pi(2);
for i = 1:ntrials
    c=rand(1);
    noh=0;
    if c<=pic
        cp=p;
    else
        cp=q;
    end        
    for j = 1:ntosses
        t=rand(1);
        if(t<=cp)
           noh=noh+1; 
        end
    end
    headcount=[headcount,noh];
end
%% Q2 EM Algorithm
est_q2=zeros(1000,3);
est_q2(1,:)=[0.5,0.5,0.45];
diff=100;
iterations=1;
while(diff>=1e-6)
    % E-Step
    sp=0;c1=[0,0];c2=[0,0];
    for i = 1:length(headcount)
        p1 = (((1-est_q2(iterations,1))/est_q2(iterations,1)) * ((est_q2(iterations,3)/est_q2(iterations,2))^headcount(i)) * (((1-est_q2(iterations,3))/(1-est_q2(iterations,2)))^(ntosses-headcount(i))) + 1); 
        p1=1/p1;
        sp=sp+p1;
        c1(1)=c1(1)+headcount(i)*p1;
        c1(2)=c1(2)+ntosses*p1;
        c2(1)=c2(1)+headcount(i)*(1-p1);
        c2(2)=c2(2)+ntosses*(1-p1);
    end

    % M-Step
    new=[0,0,0];
    new(1)=sp/ntrials;
    new(2)=c1(1)/c1(2);
    new(3)=c2(1)/c2(2);
    
    diff=max(abs(est_q2(iterations,:)-new));
    iterations=iterations+1;
    est_q2(iterations,:)=new;
end
%% PLOTS 
stem([1:1:iterations-1],est_q2(1:iterations-1,1));
xlabel('Number of iterations');
ylabel('Estimate at kth iteration');
title('Learning curve of estimate of pi');
%%
stem([1:1:iterations-1],est_q2(1:iterations-1,2));
xlabel('Number of iterations');
ylabel('Estimate at kth iteration');
title('Learning curve of estimate of p');
%%
stem([1:1:iterations-1],est_q2(1:iterations-1,3));
xlabel('Number of iterations');
ylabel('Estimate at kth iteration');
title('Learning curve of estimate of q');

%% Q3 Generating Observations
headcount=[];
ntosses=m(2);
ntrials=n(3);

a=2;
b=3;

pic=betarnd(a,b);
for i = 1:ntrials
    c=rand(1);
    noh=0;
    if c<=pic
        cp=p;
    else
        cp=q;
    end        
    for j = 1:ntosses
        t=rand(1);
        if(t<=cp)
           noh=noh+1; 
        end
    end
    headcount=[headcount,noh];
end
%% Q3 EM Algorithm
est_q3=zeros(1000,3);
iterations=1;
est_q3(iterations,:)=[0.5,0.45,0.5];
diff=100;
while(diff>=1e-6)
    % E-Step
    sp=0;c1=[0,0];c2=[0,0];
    for i = 1:length(headcount)
        p1 = (((1-est_q3(iterations,1))/est_q3(iterations,1)) * ((est_q3(iterations,3)/est_q3(iterations,2))^headcount(i)) * (((1-est_q3(iterations,3))/(1-est_q3(iterations,2)))^(ntosses-headcount(i))) + 1); 
        p1=1/p1;
        sp=sp+p1;
        c1(1)=c1(1)+headcount(i)*p1;
        c1(2)=c1(2)+ntosses*p1;
        c2(1)=c2(1)+headcount(i)*(1-p1);
        c2(2)=c2(2)+ntosses*(1-p1);
    end

    % M-Step
    new=[0,0,0];
    new(1)=sp/(ntrials+a+b-1) + (a-1)/(ntrials+a+b-1);
    new(2)=c1(1)/c1(2);
    new(3)=c2(1)/c2(2);
    
    diff=max(abs(est_q3(iterations,:)-new));
    iterations=iterations+1;
    est_q3(iterations,:)=new;
end

%% PLOTS 
stem([1:1:iterations-1],est_q3(1:iterations-1,1));
xlabel('Number of iterations');
ylabel('Estimate at kth iteration');
title('Learning curve of estimate of pi');
%%
stem([1:1:iterations-1],est_q3(1:iterations-1,2));
xlabel('Number of iterations');
ylabel('Estimate at kth iteration');
title('Learning curve of estimate of p');
%%
stem([1:1:iterations-1],est_q3(1:iterations-1,3));
xlabel('Number of iterations');
ylabel('Estimate at kth iteration');
title('Learning curve of estimate of q');


