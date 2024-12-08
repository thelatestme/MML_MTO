clear all
clc

rand('state', sum(100*clock));
alg = 1;
% func_num_single = 9;
func_num = 1;
func_end = 10;
run_times= 20;
temp = cell(1,func_end);
temp2 = cell(1,func_end);
%parpool('local',10);
tic
for index = func_num : func_end   % problem index
    [Tasks] = benchmark(index);
    maxfes = 200000;
    parfor run=1 : run_times % running times
        disp([index,run]);                    
           [population,sequence] = MM_MTO6_6_rand(Tasks, maxfes);
        temp_seq1(run,:) = sequence(:,1)';%3
        temp_seq2(run,:) = sequence(:,2)';%3
    end
    temp{alg}{index,1} = temp_seq1;%3
    temp{alg}{index,2} = temp_seq2;%3
end
time = toc;

for index = func_num : func_end   % problem index
    
    mean_data{alg}{index,1} = mean(temp{alg}{index,1},1);%3
	mean_data{alg}{index,2} = mean(temp{alg}{index,2},1);%3
    sd_data{alg}{index,1} = std(temp{alg}{index,1},1);%3
    sd_data{alg}{index,2} = std(temp{alg}{index,2},1);%3

end
for index = func_num : func_end
    curve{alg}(:,index*2-1) = mean_data{alg}{index,1};
    curve{alg}(:,index*2) = mean_data{alg}{index,2};
    mean_data1(index,1) = mean_data{alg}{index,1}(1,end);
    mean_data1(index,2) = mean_data{alg}{index,2}(1,end);
end
for alg_num=alg:alg
    data = cell(1,func_end);
    for alg=alg_num:alg_num
        for index=func_num:func_end
            for i=1:run_times
                data{index}(i,1) = temp{alg}{index,1}(i,end);
                data{index}(i,2) = temp{alg}{index,2}(i,end);
            end
        end
    end
    all_task_data = [];
    for i=1:func_end
        all_task_data = [all_task_data,data{i}];
    end
end