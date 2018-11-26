% DESCRIPTION:
%    This is the code to reproduce the result in the Fig.3

clear;
clc;
N=64; % Node Number
totalnum=N*(N-1)/2;
binary = 1;
uniform = 1;
sigma = 3*1/N;
outlier_ratio=0;
greedy_list = greedysampling(N);
for time=1:1000
    bbbbbbbbbbbbbbbbbbbbbb=time
    s=(0:N-1)'/N-1;
    s_bar = s-mean(s);
    s_bar = s_bar(randperm(N));
    index = [floor((totalnum*log(N))/N:(totalnum/100):totalnum),totalnum];

    data=[];
    for i=1:N-1
        for j=i+1:N
            data=[data
                [j i]];
        end
    end
    for i=1:totalnum
        if s_bar(data(i,1))< s_bar(data(i,2))
            temp=data(i,1);
            data(i,1)=data(i,2);
            data(i,2)=temp;
        end
    end
    i = [1:totalnum,1:totalnum];
    j = [data(:,1);data(:,2)];
    k = [ones(1,totalnum),-ones(1,totalnum)];
    d0 = full(sparse(i,j,k,totalnum,N));
    if binary==1
        if uniform ==1
            s = lsqr(d0'*d0,d0'*ones(totalnum,1)*(1-2*outlier_ratio));
        else
            s = lsqr(d0'*d0,d0'*ones(totalnum,1)*asin(1-2*outlier_ratio));
        end
        s_bar = s-mean(s);
    end


    alist=[];
    zlist=zeros(totalnum,1);

    random_sample=[]; %random sampling with replacement
    random_sample2=[];%random sampling without replacement
    active_sample=[]; %greedy sampling
    random_score=zeros(N,1);
    random_without_score=zeros(N,1);
    active_score=zeros(N,1);

    b=randperm(totalnum);
    for sample_i=1:totalnum

        if s_bar(greedy_list(sample_i,1)) < s_bar(greedy_list(sample_i,2))
            temp_1=greedy_list(sample_i,1);
            greedy_list(sample_i,1)=greedy_list(sample_i,2);
            greedy_list(sample_i,2)=temp_1;
        end
        p_re=rand(1);
        if p_re<outlier_ratio
            active_sample=[active_sample
                greedy_list(sample_i,2:-1:1)];
        else
            active_sample=[active_sample
                greedy_list(sample_i,:)];
        end

        a=ceil(totalnum*rand(1));
        alist=[alist,a];
        p_re=rand(1);
        if p_re<outlier_ratio
            random_sample=[random_sample
                data(a,2:-1:1)];
            zlist(sample_i)=1;
        else
            random_sample=[random_sample
                data(a,:)];
        end

        p_re=rand(1);
        if p_re<outlier_ratio
            random_sample2=[random_sample2
                data(b(sample_i),2:-1:1)];
        else
            random_sample2=[random_sample2
                data(b(sample_i),:)];
        end
    end

    i = [1:totalnum,1:totalnum];
    j = [random_sample(:,1);random_sample(:,2)];
    k = [ones(1,totalnum),-ones(1,totalnum)];
    d1 = full(sparse(i,j,k,totalnum,N));

    j = [random_sample2(:,1);random_sample2(:,2)];
    d2 = full(sparse(i,j,k,totalnum,N));

    j = [active_sample(:,1);active_sample(:,2)];
    d3= full(sparse(i,j,k,totalnum,N));
    epsilon1 = randn(totalnum,1);
    epsilon2 = randn(totalnum,1);
    epsilon3 = randn(totalnum,1);
    for i = 1:length(index)
        sample_i = index(i);
        if binary ==1
            if uniform ==1
                random_score = lsqr(d1(1:sample_i,:)'*d1(1:sample_i,:),d1(1:sample_i,:)'*ones(sample_i,1));
                random_without_score = lsqr(d2(1:sample_i,:)'*d2(1:sample_i,:),d2(1:sample_i,:)'*ones(sample_i,1));
                active_score = lsqr(d3(1:sample_i,:)'*d3(1:sample_i,:),d3(1:sample_i,:)'*ones(sample_i,1));
            else
                w = zeros(totalnum,1);
                z = w;
                y = w;
                for j = 1:sample_i
                    a = alist(j);
                    w(a) = w(a) + 1;
                    z(a) = z(a) + zlist(j);
                end
                y(w~=0) = asin((w(w~=0)-2*z(w~=0))./w(w~=0));
                random_score = lsqr(d0'*diag(w)*d0,d0'*diag(w)*y);
                random_without_score = lsqr(d2(1:sample_i,:)'*d2(1:sample_i,:),d2(1:sample_i,:)'*asin(ones(sample_i,1)));
                active_score = lsqr(d3(1:sample_i,:)'*d3(1:sample_i,:),d3(1:sample_i,:)'*asin(ones(sample_i,1)));
            end
        else
            random_score = s_bar +  sigma*lsqr(d1(1:sample_i,:)'*d1(1:sample_i,:),d1(1:sample_i,:)'* epsilon1(1:sample_i));
            random_without_score = s_bar +  sigma*lsqr(d2(1:sample_i,:)'*d2(1:sample_i,:),d2(1:sample_i,:)'* epsilon2(1:sample_i));
            active_score = s_bar +  sigma*lsqr(d3(1:sample_i,:)'*d3(1:sample_i,:),d3(1:sample_i,:)'* epsilon3(1:sample_i));
        end
        random_score = random_score-mean(random_score);
        random_without_score=random_without_score-mean(random_without_score);
        active_score = active_score-mean(active_score);

        % L2 distance calculation
        L2_random(i)=sqrt(sum((random_score-s_bar).^2));
        L2_random2(i)=sqrt(sum((random_without_score-s_bar).^2));
        L2_active(i)=sqrt(sum((active_score-s_bar).^2));
    end
    
    random_total(time,:)=L2_random;
    random2_total(time,:)=L2_random2;
    active_total(time,:)=L2_active;
end


for number=1:length(index)

    sd_greedy(number) = std(active_total(:,number));
    mean_greedy(number)=mean(active_total(:,number));

    sd_random(number) = std(random_total(:,number));
    mean_random(number)=mean(random_total(:,number));

    sd_random2(number) = std(random2_total(:,number));
    mean_random2(number)=mean(random2_total(:,number));
end

ratio_temp= index;
ratio_temp=ratio_temp*N/(totalnum*log(N));

 subplot(2,1,1);
plot(ratio_temp', mean_random','b','LineWidth',2 );hold on;
plot(ratio_temp',mean_random2','r','LineWidth',2 );hold on;
plot(ratio_temp',mean_greedy','g','LineWidth',2 );hold on;
legend('Random sampling with replacement','Random sampling without replacement','Greedy sampling','fontsize',24);
xlabel('p_0'); ylabel('L_2 distance');

subplot(2,1,2);
plot(ratio_temp', sd_random','b','LineWidth',2 );hold on;
plot(ratio_temp',sd_random2','r','LineWidth',2 );hold on;
plot(ratio_temp',sd_greedy','g','LineWidth',2 );hold on;
legend('Random sampling with replacement','Random sampling without replacement','Greedy sampling','fontsize',24);
xlabel('p_0'); ylabel('Standard deviation');




