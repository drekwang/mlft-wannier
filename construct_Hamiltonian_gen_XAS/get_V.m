function [V_ijkl] = get_V(i,j,k,l,V)
    if min(i,j) ~= min(k,l) || max(i,j) ~= max(k,l)
        V_ijkl = 0;
    else
        [~,idx] = ismember(V(:,1:4), [i j k l], 'rows');
        V_ijkl = V(find(idx,1),5);
    end
end