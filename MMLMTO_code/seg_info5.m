classdef seg_info5 < handle
    properties
        seg_l;%片段长度
        task;
        ind_num;
        seg_num;%原种群片段片段开始位置
        x;
        cp_task;
        cp_index;%片段所在的个体(或片段集所在的下标)
        cp_seg_index;%片段所在个体维度开始位置
        dist_min;
    end    
    methods
        function obj = seg_info5(task, ind_num, seg_num, seg_l, x)
            obj.seg_l = seg_l;
            obj.task = task;
            obj.ind_num = ind_num;
            obj.seg_num = seg_num;
            obj.x = x;
            obj.cp_task = -1;
            obj.cp_index = -1;
            obj.cp_seg_index = -1;
            obj.dist_min = inf;
        end
    end
end