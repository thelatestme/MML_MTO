function [objective] = fnceval(Task,rnvec)
    d = Task.dim;
    nvars = rnvec(1:d);
    
    minrange = Task.Lb(1:d);
    maxrange = Task.Ub(1:d);
    y=maxrange-minrange;
    vars = y.*nvars + minrange;
    x=vars;       %decoding

%     x = nvars;    %not decoding
    
    objective = Task.fnc(x);
end