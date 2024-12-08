function [population, sequence]= MM_MTO6_6_rand(Tasks, maxfes)%迁移后代全部来自交叉 %1-10，2-11，3-12
n = 75;                                                %只计算与源任务种群的最小距离
FE = 0;                                                   %每一代de之后再概率进行迁移
                                                    %不补零
seg = 0;   %片段长度 10
ntask = length(Tasks);
dims = zeros(1, ntask);
for i = 1:ntask
    dims(i) = Tasks(i).dim;
end

% result record
sequence = [];
save_gap = maxfes/10;
cur_gap = 0;

%% initilization
pop_list = {};
fit_list = {};
children_pool = cell(1,ntask);
children = cell(1,ntask);
children_fit = cell(1,ntask);
for i = 1:ntask
    %         pop = rand(n, dims(i));
    pop = (Tasks(i).Ub(1)-Tasks(i).Lb(1)).*rand(n, dims(i)) + Tasks(i).Lb(1);
    pop_list{i} = pop;
    fit = zeros(1, n);
    for j = 1:n
        fit(j)=fnceval_m(Tasks(i),pop(j, :));
        %         FE = FE + 1;
    end
    fit_list{i} = fit;
    children_pool{i}= ones(2*n,dims(i));
    children{i}= ones(n,dims(i));
    children_fit{i} = inf*ones(1,n);
end
p_kt = 0.3*ones(1,ntask);
%%
for t = 1:ntask
    pre_best_fit{t} = min(fit_list{t});
    stay_cnt{t} = 0;
     pc{t} = 0.5;
end
mtmax=0;
mtmin=1;
gen = 1;
if FE >= cur_gap
    data = [];
    for t = 1:ntask
        data = [data, min(fit_list{t})];
    end
    sequence = [sequence; data];
    disp(['fes:', num2str(FE), ' gen:', num2str(gen), ' fit1: ', num2str(data(1)), ' fit2: ', num2str(data(2))]);
    cur_gap = cur_gap + save_gap;
end
while FE < maxfes

    for t = 1:ntask
        if pre_best_fit{t} == min(fit_list{t})
            stay_cnt{t} = stay_cnt{t}+1;
            if stay_cnt{t} >= 2
                stay_cnt{t}= 0;
                pc{t} = rand()*(mtmax-mtmin)+mtmin;
            end
        else
            pre_best_fit{t} = min(fit_list{t});
            stay_cnt{t}= 0;
        end
    end

    gen = gen + 1;
    for t = 1:ntask
        cur_pop = pop_list{t};
        cur_fit = fit_list{t};
           
            for i = 1:n
                r1 = randomInt(1, n);
                r2 = randomInt(1, n);
                r3 = randomInt(1, n);
                while r1 == r2 || r2 == r3 || r3 == r1 || r1 == i || r2 == i || r3 == i
                    r1 = randomInt(1, n);
                    r2 = randomInt(1, n);
                    r3 = randomInt(1, n);
                end
                v = cur_pop(r1, :) + 0.5*(cur_pop(r2, :) - cur_pop(r3, :));
                u = cur_pop(i, :);
                jrand = randomInt(1, dims(t));
                for j = 1:dims(t)
                    if rand() < 0.9 ||  j == jrand
                        u(j) = v(j);
                    end
                end
                %             u = max(0, min(1, u));
                u = max(Tasks(t).Lb, min(Tasks(t).Ub, u));
                ufit = fnceval_m(Tasks(t),u);
                FE = FE + 1;
                if FE >= cur_gap && FE<maxfes
                    data = [];
                    for ti = 1:ntask
                        data = [data, min(fit_list{ti})];
                    end
                    sequence = [sequence; data];
                    disp(['fes:', num2str(FE), ' gen:', num2str(gen), ' fit1: ', num2str(data(1)), ' fit2: ', num2str(data(2))]);
                    cur_gap = cur_gap + save_gap;
                end
                if ufit <= cur_fit(i)
                    cur_pop(i, :) = u;
                    cur_fit(i) = ufit;
                end
                if FE>=maxfes
                    break;
                end
            end

        ktn = 5;%只对前knt个个体进行迁移     
        if rand()< p_kt                           
            st = randi(ntask);                     
            while st == t          
                st = randi(ntask);
            end
            %排序
            [fit_list{t},sort_index] = sort(fit_list{t},'ascend');
            pop_list{t} = pop_list{t}(sort_index,:);
            [fit_list{st},sort_index] = sort(fit_list{st},'ascend');
            pop_list{st} = pop_list{st}(sort_index,:);

            info_num = 0;
            cur_pop = pop_list{t};
            %分割当前种群为片段
            for class=0:2
                for i=1:ktn
                    seg = floor(rand()*(16-5)+5);%floor(rand()*(7-2)+2);
                    dim = dims(t);
                    if mod(dims(t),seg)~=0
                        dim = ceil(dims(t)/seg)*seg;
                    end
                    for p = 1:seg:dim
                        info_num = info_num+1;
                        try
                            x = cur_pop(class*25+i,p:p+seg-1);
                            temp_seg = seg;
                        catch
                            gap = dim-dims(t);%维度不够,不处理，直接为小片段
                            x = cur_pop(class*25+i,p:end);
                            temp_seg = seg-gap;
                        end
                        seg_list{t}(info_num) = seg_info5(t,class*25+i,p,temp_seg,x);
                    end
                end
            end
          
            for ii=1:length(seg_list{t})
                
                temp_seg = seg_list{t}(ii).seg_l;
                ni = seg_list{t}(ii).ind_num;
                class = floor(ni/25);
                for i = 1 : 15
                    seg_num = dims(st)-temp_seg+1;
                    for j=1:seg_num
                        st_x = pop_list{st}(i, j:j+temp_seg-1);
                        dist = norm(seg_list{t}(ii).x - st_x);
                        if dist < seg_list{t}(ii).dist_min
                            seg_list{t}(ii).cp_index = i;
                            seg_list{t}(ii).cp_seg_index = j;
                            seg_list{t}(ii).cp_task = st;
                            seg_list{t}(ii).dist_min = dist;
                        end
                    end
                end
                %欧式距离最小的两个片段进行交叉
                cp_i = seg_list{t}(ii).cp_index;
                cp_j = seg_list{t}(ii).cp_seg_index;
                cp_t = seg_list{t}(ii).cp_task;
                cp_x = pop_list{cp_t}(cp_i, cp_j:cp_j+temp_seg-1);
                dist_min = seg_list{t}(ii).dist_min;
                x = seg_list{t}(ii).x;
                %迁移交流
                sel_arr = [0.3,0.2,pc{t}];
                sel = RouletteWheelSelection(sel_arr); 					
                switch sel
                    case 1
                        cr = randomDouble(0.1, 0.9);        %方法：二项式交叉
                        drand = randomInt(1, temp_seg);
                        for d = 1:temp_seg
                            if rand() < cr || d == drand
                                x(d) = cp_x(d);
                            end
                        end
                    case 2
                        ud = rand(1,temp_seg);                 %方法：SBX  (偏向cp_x
                        cf = zeros(1,temp_seg);
                        mu1 = 2;
                        cf(ud<=0.5)=(2*ud(ud<=0.5)).^(1/(mu1+1));
                        cf(ud>0.5)=(2*(1-ud(ud>0.5))).^(-1/(mu1+1));
                        x = SBX_inside(cp_x,x,cf);
                    case 3
                        sigma = dist_min/temp_seg;
                        rvec = normrnd(0,sigma*0.01,[1,temp_seg]);
                        x = x+rvec;
                end
                seg_t = seg_list{t}(ii).task;
                seg_i = seg_list{t}(ii).ind_num;
                seg_n = seg_list{t}(ii).seg_num;
                %                 x = max(0, min(1, x));
                x = max(Tasks(seg_t).Lb(1), min(Tasks(seg_t).Ub(1), x));%过界处理 
                children_pool{seg_t}(seg_i,seg_n:seg_n+temp_seg-1) = x;
            end
            %评估后代适应值
            dim = size(children_pool{t},2);
            if  dim > dims(t)
                children_pool{t}(:,dims(t)+1:dim) = [];
            end
            children{t} = children_pool{t};
            for class=0:2
                for i=1:ktn
                    children_fit{t}(class*25+i) = fnceval_m(Tasks(t),children{t}(class*25+i,:));
                    FE = FE+1;
                    if FE >= cur_gap && FE<maxfes
                        data = [];
                        for ti = 1:ntask
                            data = [data, min(fit_list{ti})];
                        end
                        sequence = [sequence; data];
                        disp(['fes:', num2str(FE), ' gen:', num2str(gen), ' fit1: ', num2str(data(1)), ' fit2: ', num2str(data(2))]);
                        cur_gap = cur_gap + save_gap;
                    end
                end
            end

            for class=0:2           %%25为同一等级的个体数 25+ktn
                temp_fit_p = fit_list{t}(class*25+1 : (class+1)*25);       
                temp_fit_c = children_fit{t}(class*25+1 : class*25+ktn);
                temp_p = pop_list{t}(class*25+1 : (class+1)*25 , :);
                temp_c = children{t}(class*25+1 : class*25+ktn , :);
                class_pop = [temp_p;temp_c];
                class_fit = [temp_fit_p,temp_fit_c];
                [class_fit,sort_index] = sort(class_fit,'ascend');
                cur_pop(class*25+1 : (class+1)*25 , :) = class_pop(sort_index(1:25),:);
                cur_fit(class*25+1 : (class+1)*25) = class_fit(1:25);
            end

        end
        %更新种群
        pop_list{t} = cur_pop;
        fit_list{t} = cur_fit;
    end

end

%%
data = [];
for t = 1:ntask
    data = [data, min(fit_list{t})];
end
sequence = [sequence; data];
disp(['fes:', num2str(FE), ' gen:', num2str(gen), ' fit1: ', num2str(data(1)), ' fit2: ', num2str(data(2))]);

for t = 1:ntask
    [fit_list{t},sort_index] = sort(fit_list{t},'ascend');
    pop_list{t} = pop_list{t}(sort_index,:);
end
population = pop_list;
end

function object=SBX_inside(p1,p2,cf)
object=0.5*((1+cf).*p1 + (1-cf).*p2);
% object(object>1)=1;
% object(object<0)=0;
end