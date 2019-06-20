function [r_range, R] = gen_radialfunctions(...
    xsffile, npts, lengths, origin,atomcenter)
[x_mesh, y_mesh, z_mesh, cartcoords] = gen_cartcoords(npts,lengths,origin,atomcenter);
[values, values_3D] = gen_values(xsffile, npts);
[sphcoords] = gen_sphcoords(cartcoords);
[sphcoords_ordered,values_ordered]=order_sphscalarfield(sphcoords, values);
%Limiting makes the code run way faster because the interpolation step is 
%super slow. For SrTiO3, r<2.3 A at 40,000 points, r<3 at 88603 points, r<4
%at 209971, r<5 at 410018. 
%Test: No real difference in F_k between r<3 and r<4 for SrTiO3_00001
%(dz^2)
index_cutoff = 209971;
sphcoords_ordered = sphcoords_ordered(1:index_cutoff,:);
values_ordered = values_ordered(1:index_cutoff,:);
[r_range, R] = integrate_angular(sphcoords_ordered,values_ordered);
end



