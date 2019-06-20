function [qn_index] = get_qn_index(n,l,m,spin,atom,qn)
    [~,idx] = ismember(qn(:,1:5),...
            [n l m spin atom], 'rows');
    %There should only be one match for '[l1 m1 l2 m2 order]'
    qn_index = find(idx,1);
end