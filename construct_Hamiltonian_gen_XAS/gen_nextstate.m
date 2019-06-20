function [nextstate] = gen_nextstate(state)
    tempstate = state;
    
    [~,index] = max(state);
    if index > 1
        tempstate(1,index) = 0;
        tempstate(1,index-1) = 1;
    else
        [minvalue,minindex] = min(state);
        numatend = minindex - 1;
        for i = 1:numatend
            tempstate(1,i) = 0;
        end
        [~,index] = max(tempstate);
        tempstate(1,index) = 0;
        for i = index-1-numatend:index-1
            tempstate(1,i) = 1;
        end
    end 
    nextstate = tempstate;
end