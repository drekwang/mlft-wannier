function [Rk_ijkl] = get_R(ni,li,nj,lj,nk,lk,nl,ll,atomk,atoml,order,Rk)
    if atoml == 0 || atomk == 0
        Rk_ijkl = 0;
    else
        [~,idx] = ismember(Rk(:,1:11),...
            [ni li nj lj nk lk nl ll atomk atoml order], 'rows');
        Rk_ijkl = Rk(find(idx,1),12);
    end
end