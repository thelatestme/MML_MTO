classdef seg_info5 < handle
    properties
        seg_l;%Ƭ�γ���
        task;
        ind_num;
        seg_num;%ԭ��ȺƬ��Ƭ�ο�ʼλ��
        x;
        cp_task;
        cp_index;%Ƭ�����ڵĸ���(��Ƭ�μ����ڵ��±�)
        cp_seg_index;%Ƭ�����ڸ���ά�ȿ�ʼλ��
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