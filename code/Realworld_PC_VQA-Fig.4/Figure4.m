% DESCRIPTION:
%   This is the code to reproduce the result in the Fig.4.
 

clc;
clear;
warning off;

data_ind=[];
N=16;
totalnum=N*(N-1)/2;
index_initial=floor((120*log(N))/N);
greedy_list = greedysampling(N);
ref_num=10;
for i=1:N-1
    for j=i+1:N
        data_ind=[data_ind
            [j i]];
    end
end

for ref=1:ref_num; %reference index, there are 10 references in PC-VQA dataset
    totalround=32; % 32 round complete paired comparison in PC-VQA dataset
    str=strcat('.\PC-VQA dataset\','data',num2str(ref),'.mat');
    load(str); % load reference  
    %  Input matrix with paired comparison, the matrix has
    %  2 columns, and for each row, the rank of the first column
    %  is higher than the second column for this comparison.
    ground_truth_score=Hodgerank(data_ref);
    s_bar=ground_truth_score-mean(ground_truth_score);


    for time=1:100
        greedy_perm=randperm(N);

        for ii=1:totalnum
            greedy_list(ii,1)=greedy_perm(greedy_list(ii,1));
            greedy_list(ii,2)=greedy_perm(greedy_list(ii,2));

        end
        a=randperm(totalnum);

        active_sample=[];
        random_sample=[];
        random_sample2=[];

        L2_active=zeros(1,totalnum);
        L2_random=zeros(1,totalnum);
        L2_random2=zeros(1,totalnum);


        for sample_i=1:totalnum
            active_score=zeros(N,1);
            random_score=zeros(N,1);
            random_without_score=zeros(N,1);

            b=randperm(totalnum);

            for i=1:totalnum
                round=randperm(totalround);
                round_select1=round(1);
                j=(round_select1-1)*totalnum+1:round_select1*totalnum;
                data=data_ref(j,:);
                if data(i,1)==data_ind(a(sample_i),1) && data(i,2)==data_ind(a(sample_i),2) || data(i,1)==data_ind(a(sample_i),2)&& data(i,2)==data_ind(a(sample_i),1)
                    random_sample=[random_sample
                        data(i,:)];
                end
                round=randperm(totalround);
                round_select2=round(1);
                j=(round_select2-1)*totalnum+1:round_select2*totalnum;
                data=data_ref(j,:);
                if data(i,1)==data_ind(b(sample_i),1) && data(i,2)==data_ind(b(sample_i),2) || data(i,1)==data_ind(b(sample_i),2)&& data(i,2)==data_ind(b(sample_i),1)
                    random_sample2=[random_sample2
                        data(i,:)];
                end
                round=randperm(totalround);
                round_select3=round(1);
                j=(round_select3-1)*totalnum+1:round_select3*totalnum;
                data=data_ref(j,:);
                if data(i,1)==greedy_list(sample_i,1) && data(i,2)==greedy_list(sample_i,2) || data(i,1)==greedy_list(sample_i,2) && data(i,2)==greedy_list(sample_i,1)
                    active_sample=[active_sample
                        data(i,:)];
                end
            end

            if sample_i>=floor((totalnum*log(N))/N)
                active_score=Hodgerank(active_sample);
                active_score=active_score-mean(active_score);
            end


            if sample_i>=floor((totalnum*log(N))/N)
                random_without_score=Hodgerank(random_sample);
                random_without_score=random_without_score-mean(random_without_score);
            end

            if sample_i>=floor((totalnum*log(N))/N)
                random_score=Hodgerank(random_sample2);
                random_score=random_score-mean(random_score);
            end

            if length(random_without_score)==length(s_bar)&& length(s_bar)==length(active_score)&& length(random_score)==length(s_bar)
                L2_active(sample_i)=sqrt(sum((active_score-s_bar).^2)); % l2 distance
                L2_random(sample_i)=sqrt(sum((random_without_score-s_bar).^2));
                L2_random2(sample_i)=sqrt(sum((random_score-s_bar).^2));
            end

            L2_random(1:index_initial)=NaN;
            L2_active(1:index_initial)=NaN;
            L2_random2(1:index_initial)=NaN;
        end

        active_total(time,:)=L2_active;
        random_total(time,:)=L2_random;
        random2_total(time,:)=L2_random2;

    end


    for number=1:totalnum

        mean_greedy(number)=mean(active_total(:,number));
        mean_random(number)=mean(random_total(:,number));
        mean_random2(number)=mean(random2_total(:,number));

    end


    ratio_temp= 1:number;
    ratio_temp=ratio_temp*N/(totalnum*log(N));
    cutnum = sum(ratio_temp<1);
%     plot(ratio_temp',mean_random','r','LineWidth',2 );hold on;
%     plot(ratio_temp',mean_random2','b','LineWidth',2 );hold on;
%     plot(ratio_temp',mean_greedy','g','LineWidth',2 );hold on;
% 
%     legend('Random sampling without replacement','Random sampling with replacement','Greedy sampling','fontsize',24);
%     xlabel('p_0'); ylabel('L_2 distance');
    str_result=strcat('ref',num2str(ref),'.mat');
    save(str_result);
end

% show the results

for ref=1:ref_num
    str=strcat('ref',num2str(ref),'.mat')
    load(str);
    subplot(3,4,ref);
    for number=1:totalnum
        mean_greedy(number)=mean(active_total(:,number));
        mean_random(number)=mean(random_total(:,number));
        mean_random2(number)=mean(random2_total(:,number));
    end
    ratio_temp= 1:number;
    ratio_temp=ratio_temp*N/(totalnum*log(N));
    cutnum = sum(ratio_temp<1);
    plot(ratio_temp',mean_random','r','LineWidth',2 );hold on;
    plot(ratio_temp',mean_random2','b','LineWidth',2 );hold on;
    plot(ratio_temp',mean_greedy','g','LineWidth',2 );hold on;
    if ref==ref_num
    legend('Random sampling without replacement','Random sampling with replacement','Greedy sampling','fontsize',24);
    end
    xlabel('p_0'); ylabel('L_2 distance');
end
