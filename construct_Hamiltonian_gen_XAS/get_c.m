function [c] = get_c(l1,m1,l2,m2,order,c)
    [~,idx] = ismember(c(:,1:5),...
            [l1 m1 l2 m2 order], 'rows');
    %There should only be one match for '[l1 m1 l2 m2 order]'
    c = c(find(idx,1),6);
end