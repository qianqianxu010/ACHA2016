% DESCRIPTION:
%    This is the code to reproduce the result in the Fig.2

clear;clc;

N_list = [32,128,512,2048];
for ii = 1:4
    subplot(2,2,ii);
    data=[];
    N=N_list(ii);
    
    totalpair=N*(N-1)/2;
    for i=1:N-1  %% the complete data
        data=[data;[i*ones(N-i,1),(i+1:N)']];
    end
    
    p0_max = N/log(N);  %% plot points
    p0 = 1:(N/log(N)-1)/30:p0_max;
    ind = floor(p0/2*(N-1)*log(N));
    ind(length(p0))=totalpair;
    iter=20;
    
    lamda_with=zeros(iter,length(p0));% %  lambda, with replacement
    lamda_without=zeros(iter,length(p0));%  lambda, without replacement
    dg_with = zeros(iter,length(p0));   %%  degree, with replacement
    dg_without = zeros(iter,length(p0));   %%  degree,  without replacement
    for time=1:iter
        bb=randperm(totalpair);   %% random sampling without replacement
        a=ceil(totalpair*rand(max(ind),1));  %%%% random sampling with replacement
        
        i = [1:max(ind),1:max(ind)];
        j = [data(a,1);data(a,2)];
        k = [ones(1,max(ind)),-ones(1,max(ind))];
        d1 = sparse(i,j,k,max(ind),N);
        j = [data(bb(1:max(ind)),1);data(bb(1:max(ind)),2)];
        d2 = sparse(i,j,k,max(ind),N);
        for k = 1:length(p0)
            L = d1(1:ind(k),:)'*d1(1:ind(k),:);
            dg1(k) = min(sum(abs(d1(1:ind(k),:))))/log(N)/p0(k);
            [V,D] = eigs(L,[],2,'SA');
            lamda1(k)=D(2,2)/log(N)/p0(k);
            
            L = d2(1:ind(k),:)'*d2(1:ind(k),:);
            dg2(k) = min(sum(abs(d2(1:ind(k),:))))/log(N)/p0(k);
            [V,D] = eigs(L,[],2,'SA');
            lamda2(k)=D(2,2)/log(N)/p0(k);
        end
        lamda_with(time,:)=lamda1;
        lamda_without(time,:)=lamda2;
        dg_with(time,:)=dg1;
        dg_without(time,:)=dg2;
    end
    for number=1:length(p0)
        %quan_lambda_with=quantile(lamda_with(:,number),[0.25,0.5,0.75]);
        %quan_lambda_without=quantile(lamda_without(:,number),[0.25,0.5,0.75]);
        mean_lambda_with(number)=mean(lamda_with(:,number));
        mean_lambda_without(number)=mean(lamda_without(:,number));
        mean_dg_with(number)=mean(dg_with(:,number));
        mean_dg_without(number)=mean(dg_without(:,number));
    end
    
    plot(p0,mean_lambda_with,'b',p0,mean_lambda_without,'r','LineWidth',2);
    hold on;
    plot(p0,mean_dg_with,'b.',p0,mean_dg_without,'r.','LineWidth',2)
    ezplot('x-1-y*x*(1-log(y))',[1,max(p0),0,1]);
    plot(p0, 1-sqrt(2./p0.*(1-p0*log(N)/N)),'m');
    
    title('')
    ylim([0,1]);
    xlabel('p_0')
    ylabel('\lambda_2/np or d_{min}/np')
    
    legend('Random sampling with replacement','Random sampling without replacement',...
        'Minimal degree of with replacement', 'Minimal degree of without replacement', ...
        'a(p_0)','a_2(p_0)','fontsize',20);
end