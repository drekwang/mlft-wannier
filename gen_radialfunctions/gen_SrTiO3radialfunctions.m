numfunctions = 14;
R_Wannier = zeros(100,numfunctions);
% for i = 1:5
%     filename = strcat(strcat('SrTiO3_',num2str(i)),'.xsf');
%     [r_Wannier,R]=gen_radialfunctions(filename,npts,lengths,origin,atomcenter);
%     R_Wannier(:,i) = R;
% end
atomcenter = [0 1.9525 1.9525];
for i = 6:8
    filename = strcat(strcat('SrTiO3_',num2str(i)),'.xsf');
    [r_Wannier,R]=gen_radialfunctions(filename,npts,lengths,origin,atomcenter);
    R_Wannier(:,i) = R;
end
% atomcenter = [1.9525 1.9525 0];
% for i = 9:11
%     filename = strcat(strcat('SrTiO3_',num2str(i)),'.xsf');
%     [r_Wannier,R]=gen_radialfunctions(filename,npts,lengths,origin,atomcenter);
%     R_Wannier(:,i) = R;
% end
% atomcenter = [1.9525 0 1.9525];
% for i = 12:14
%     filename = strcat(strcat('SrTiO3_',num2str(i)),'.xsf');
%     [r_Wannier,R]=gen_radialfunctions(filename,npts,lengths,origin,atomcenter);
%     R_Wannier(:,i) = R;
% end