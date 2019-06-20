function [r,R,area] = read_LDA_AE(file, n, l)
    %'Ti_LDA_AE.orbs'
    %l=0 for s, l=1 for p, l=2 for d    
    data = importdata(file,' ');
    R_cutoff = 0.01;
    index = and(any(data(:,1) == n, 2), and(any(data(:,2) == l, 2), any(abs(data(:,4))>R_cutoff,2)));
    index2 = find(index==1);
    count = 1;
    r = zeros(size(index2,1),1);
    R = zeros(size(index2,1),1);
    for i = index2(1,1):index2(end,1)
        r(count) = data(i,3);
        %What's in the file is r*R, so divide by r
        R(count) = data(i,4)/data(i,3);
        count = count+1;
    end
    %plot(r,R)
    area = trapz(r,R.^2.*r.^2);
end