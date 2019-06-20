% %Check change in theta/phi convention
% sphcoords_i = zeros(size(sphcoords_ordered,1),3);
% for i = 1:size(sphcoords_ordered,1)
%     if sphcoords_ordered(i,1)>pi()
%         sphcoords_i(i,1)=sphcoords_ordered(i,1)-2*pi();
%     else
%         sphcoords_i(i,1)=sphcoords_ordered(i,1);
%     end
%     sphcoords_i(i,2) = pi()/2 - sphcoords_ordered(i,2);
% end
% sphcoords_i(:,3) = sphcoords_ordered(:,3);
% [x y z]=sph2cart(sphcoords_i(:,1),sphcoords_i(:,2),sphcoords_i(:,3));
% figure; scatter3(x,y,z,abs(values_ordered)); xlabel('x');ylabel('y');zlabel('z')

% %Check griddata interpolation
% sphpoints = [phi_mesh(:) theta_mesh(:) r_mesh(:) values_mesh(:)]; 
% sphcoords_i = zeros(size(sphpoints,1),4);
% for i = 1:size(sphpoints,1)
%     sphcoords_i(i,1) = sphpoints(i,1);
%     sphcoords_i(i,2) = pi()/2 - sphpoints(i,2);
% end
% sphcoords_i(:,3) = sphpoints(:,3);
% sphcoords_i(:,4) = sphpoints(:,4);
% [x_sph y_sph z_sph]=sph2cart(sphcoords_i(:,1)',sphcoords_i(:,2)',sphcoords_i(:,3)'); 
% close all; 
% scatter3(x_sph',y_sph',z_sph',abs(sphcoords_i(:,4))+0.00001);
% xlabel('x');ylabel('y');zlabel('z')

% %Check mirroring
% sphcoords_i = zeros(size(mirrored,1),3);
% for i = 1:size(mirrored,1)
%     sphcoords_i(i,1) = mirrored(i,1);
%     sphcoords_i(i,2) = pi()/2 - mirrored(i,2);
% end
% sphcoords_i(:,3) = mirrored(:,3);
% [x y z]=sph2cart(sphcoords_i(:,1),sphcoords_i(:,2),sphcoords_i(:,3));
% figure; scatter3(x,y,z,abs(mirrored(:,4))); xlabel('x');ylabel('y');zlabel('z')

% %Integrate wavefunction in Cartesian coordinates to see if it's normalized
% cartpts = [linspace(origin(1)-atomcenter(1),origin(1)+lengths(1)-atomcenter(1),npts(1));...
%     linspace(origin(2)-atomcenter(2),origin(2)+lengths(2)-atomcenter(2),npts(2));...
%     linspace(origin(3)-atomcenter(3),origin(3)+lengths(3)-atomcenter(3),npts(3))]';
% x_range = cartpts(:,1); y_range = cartpts(:,2); z_range = cartpts(:,3);
% cartpts = [linspace(0,lengths(1),npts(1));...
%     linspace(0,lengths(2),npts(2));...
%     linspace(0,lengths(3),npts(3))]';
% x_range = cartpts(:,1); y_range = cartpts(:,2); z_range = cartpts(:,3);
% %int = trapz(z_range,trapz(x_range,trapz(y_range,values_3D.^2,2)))
% int1 = trapz(x_range,values_3D.^2,3);
% int2 = trapz(y_range,int1,2);
% int3 = trapz(z_range,int2,1)

% %Integrate only 1/2 of the wavefunction in Cartesian coordinates to see if
% %we've centered it correctly (z>0)
% %Apparently, the simulation cell does not have the same bounds in the
% %+-x/y/z directions: -5.9635 vs. 5.6431, so cutting numpoints in half
% %won't halve the integral
% cartcoords1 = cartcoords(1:size(cartcoords,1)/2,:);
% values1 = values(1:size(cartcoords,1)/2,1)';
% values_3D1 = zeros(108,108,54);
% count = 1;
% for k = 1:54
%     for j = 1:108
%         for i = 1:108
%             values_3D1(i,j,k) = values1(1,count);
%             count = count+1;
%         end
%     end
% end
% values1 = values1';
% x_range = linspace(origin(1)-atomcenter(1),origin(1)+lengths(1)-atomcenter(1),108); 
% y_range = linspace(origin(2)-atomcenter(2),origin(2)+lengths(2)-atomcenter(2),108); 
% z_range = linspace(0,origin(3)+lengths(3)-atomcenter(3),54);
% int = trapz(z_range,trapz(x_range,trapz(y_range,values_3D1.^2,2)))

% %Plot Wannier function using Cartesian coordinates
% scatter3(x_mesh(:),y_mesh(:),z_mesh(:),abs(values_3D(:)));
% hold on
% scatter3(0,0,0,100,'r')

% %Test code for Slater integration
% %Constants
% e = 1.6*10^-19; %elementary charge in C
% eps_0 = 8.854187*10^-12; %vacuum permittivity in F/m = s^2*C^2*m^-2*kg^-1
% %input
% k = 2;
% r1 = r_plot;
% r2 = r_plot;
% R1 = R;
% R2 = R;
% [r1_mesh r2_mesh] = meshgrid(r1,r2);
% r_points = [r1_mesh(:) r2_mesh(:)];
% values = zeros(size(r_points,1),1);
% for i = 1:size(r_points,1)
%     values(i,1) = min(r_points(i,1),r_points(i,2))^k/...
%         (max(r_points(i,1),r_points(i,2)))^(k+1)*r_points(i,1)^2*...
%         r_points(i,2)^2*R1(find(r1==r_points(i,1)))^2*R2(find(r2==r_points(i,2)))^2;
% end
% all = [r_points values];
% all(isnan(all))=0;
% values_mesh = zeros(size(r_plot,1),size(r_plot,1));
% count = 1;
% for j = 1:size(r_plot,1)
%     for i = 1:size(r_plot,1)
%         values_mesh(i,j) = all(count,3);
%         count = count+1;
%     end
% end
% %F_k in atomic units (Angstrom^-1)
% F_k_au = trapz(r2,trapz(r1,values_mesh,2))
% %Convert from atomic units to eV
% F_k = F_k_au*e^2/(4*pi()*eps_0)/e/10^-10
sum = 0;
for i = 3:52
    sum = sum+HDFT(i,i);
end