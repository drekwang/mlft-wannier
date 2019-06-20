%The loaded input file should have 'atomcenter' as 1x3 array, 'lengths'
%as 1x3 array, 'npts' as 1x3 array, and 'origin' as 1x3 array
load('NiO_Ni_input.mat')
%Need to integrate 5 3d and 1 2p
xsffile1 = 'NiO_00001.xsf'; %d_z^2
xsffile2= 'NiO_00002.xsf'; %d_xz
xsffile3 = 'NiO_00003.xsf'; %d_yz
xsffile4 = 'NiO_00004.xsf'; %d_x^2-y^2
xsffile5 = 'NiO_00005.xsf'; %d_xy

%Generate 'x/y/z_mesh' and a list of coordinates 'cartcoords' using number 
%of points in each x/y/z direction 'npts', length of the unit cell in each
%direction 'lengths', the origin of the box 'origin', and where the atom of
%interested in located 'atomcenter'
[x_mesh, y_mesh, z_mesh, cartcoords] = gen_cartcoords(...
    npts,lengths,origin,atomcenter);
fprintf('Created x/y/z mesh and cartcoords')
fprintf('\n')
%'values' is a list of values of the wavefunction in the same order as 
%'cartcoords'; 'values_3D' is a 3D matrix to plot with 'x/y/z_mesh'
[values, values_3D] = gen_values(xsffile1, npts);
fprintf('Loaded in values from xsf')
fprintf('\n')
%Convert cartesian coordinates to spherical coordinates 'sphcoords'
[sphcoords] = gen_sphcoords(cartcoords);
fprintf('Converted from Cartesian to spherical coordinates')
fprintf('\n')
%Order the spherical coordinates to make it easier to integrate
[sphcoords_ordered,values_ordered]=order_sphscalarfield(sphcoords, values);
fprintf('Ordered spherical coordinates and corresponding values')
fprintf('\n')
%Limiting the 'index_cutoff'makes the code run way faster because the 
%interpolation step is super slow. For NiO, r<3 A at 99562 points
%Test: No real difference in F_k between r<3 and r<4 for NiO_00001 (dz^2)
index_cutoff = 99562;
sphcoords_ordered = sphcoords_ordered(1:index_cutoff,:);
values_ordered = values_ordered(1:index_cutoff,:);
%Integrate the wavefunction across theta and phi to get the radial wavefn R
%as a function of 'r_range'
fprintf('Beginning angular integration of wavefunction')
fprintf('\n')

[r_range1, R1] = integrate_angular(sphcoords_ordered,values_ordered);

%Generate 'x/y/z_mesh' and a list of coordinates 'cartcoords' using number 
%of points in each x/y/z direction 'npts', length of the unit cell in each
%direction 'lengths', the origin of the box 'origin', and where the atom of
%interested in located 'atomcenter'
[x_mesh, y_mesh, z_mesh, cartcoords] = gen_cartcoords(...
    npts,lengths,origin,atomcenter);
fprintf('Created x/y/z mesh and cartcoords')
fprintf('\n')
%'values' is a list of values of the wavefunction in the same order as 
%'cartcoords'; 'values_3D' is a 3D matrix to plot with 'x/y/z_mesh'
[values, values_3D] = gen_values(xsffile2, npts);
fprintf('Loaded in values from xsf')
fprintf('\n')
%Convert cartesian coordinates to spherical coordinates 'sphcoords'
[sphcoords] = gen_sphcoords(cartcoords);
fprintf('Converted from Cartesian to spherical coordinates')
fprintf('\n')
%Order the spherical coordinates to make it easier to integrate
[sphcoords_ordered,values_ordered]=order_sphscalarfield(sphcoords, values);
fprintf('Ordered spherical coordinates and corresponding values')
fprintf('\n')
%Limiting the 'index_cutoff'makes the code run way faster because the 
%interpolation step is super slow. For NiO, r<3 A at 99562 points
%Test: No real difference in F_k between r<3 and r<4 for NiO_00001 (dz^2)
index_cutoff = 99562;
sphcoords_ordered = sphcoords_ordered(1:index_cutoff,:);
values_ordered = values_ordered(1:index_cutoff,:);
%Integrate the wavefunction across theta and phi to get the radial wavefn R
%as a function of 'r_range'
fprintf('Beginning angular integration of wavefunction')
fprintf('\n')

[r_range2, R2] = integrate_angular(sphcoords_ordered,values_ordered);

%Generate 'x/y/z_mesh' and a list of coordinates 'cartcoords' using number 
%of points in each x/y/z direction 'npts', length of the unit cell in each
%direction 'lengths', the origin of the box 'origin', and where the atom of
%interested in located 'atomcenter'
[x_mesh, y_mesh, z_mesh, cartcoords] = gen_cartcoords(...
    npts,lengths,origin,atomcenter);
fprintf('Created x/y/z mesh and cartcoords')
fprintf('\n')
%'values' is a list of values of the wavefunction in the same order as 
%'cartcoords'; 'values_3D' is a 3D matrix to plot with 'x/y/z_mesh'
[values, values_3D] = gen_values(xsffile3, npts);
fprintf('Loaded in values from xsf')
fprintf('\n')
%Convert cartesian coordinates to spherical coordinates 'sphcoords'
[sphcoords] = gen_sphcoords(cartcoords);
fprintf('Converted from Cartesian to spherical coordinates')
fprintf('\n')
%Order the spherical coordinates to make it easier to integrate
[sphcoords_ordered,values_ordered]=order_sphscalarfield(sphcoords, values);
fprintf('Ordered spherical coordinates and corresponding values')
fprintf('\n')
%Limiting the 'index_cutoff'makes the code run way faster because the 
%interpolation step is super slow. For NiO, r<3 A at 99562 points
%Test: No real difference in F_k between r<3 and r<4 for NiO_00001 (dz^2)
index_cutoff = 99562;
sphcoords_ordered = sphcoords_ordered(1:index_cutoff,:);
values_ordered = values_ordered(1:index_cutoff,:);
%Integrate the wavefunction across theta and phi to get the radial wavefn R
%as a function of 'r_range'
fprintf('Beginning angular integration of wavefunction')
fprintf('\n')

[r_range3, R3] = integrate_angular(sphcoords_ordered,values_ordered);

%Generate 'x/y/z_mesh' and a list of coordinates 'cartcoords' using number 
%of points in each x/y/z direction 'npts', length of the unit cell in each
%direction 'lengths', the origin of the box 'origin', and where the atom of
%interested in located 'atomcenter'
[x_mesh, y_mesh, z_mesh, cartcoords] = gen_cartcoords(...
    npts,lengths,origin,atomcenter);
fprintf('Created x/y/z mesh and cartcoords')
fprintf('\n')
%'values' is a list of values of the wavefunction in the same order as 
%'cartcoords'; 'values_3D' is a 3D matrix to plot with 'x/y/z_mesh'
[values, values_3D] = gen_values(xsffile4, npts);
fprintf('Loaded in values from xsf')
fprintf('\n')
%Convert cartesian coordinates to spherical coordinates 'sphcoords'
[sphcoords] = gen_sphcoords(cartcoords);
fprintf('Converted from Cartesian to spherical coordinates')
fprintf('\n')
%Order the spherical coordinates to make it easier to integrate
[sphcoords_ordered,values_ordered]=order_sphscalarfield(sphcoords, values);
fprintf('Ordered spherical coordinates and corresponding values')
fprintf('\n')
%Limiting the 'index_cutoff'makes the code run way faster because the 
%interpolation step is super slow. For NiO, r<3 A at 99562 points
%Test: No real difference in F_k between r<3 and r<4 for NiO_00001 (dz^2)
index_cutoff = 99562;
sphcoords_ordered = sphcoords_ordered(1:index_cutoff,:);
values_ordered = values_ordered(1:index_cutoff,:);
%Integrate the wavefunction across theta and phi to get the radial wavefn R
%as a function of 'r_range'
fprintf('Beginning angular integration of wavefunction')
fprintf('\n')

[r_range4, R4] = integrate_angular(sphcoords_ordered,values_ordered);

%Generate 'x/y/z_mesh' and a list of coordinates 'cartcoords' using number 
%of points in each x/y/z direction 'npts', length of the unit cell in each
%direction 'lengths', the origin of the box 'origin', and where the atom of
%interested in located 'atomcenter'
[x_mesh, y_mesh, z_mesh, cartcoords] = gen_cartcoords(...
    npts,lengths,origin,atomcenter);
fprintf('Created x/y/z mesh and cartcoords')
fprintf('\n')
%'values' is a list of values of the wavefunction in the same order as 
%'cartcoords'; 'values_3D' is a 3D matrix to plot with 'x/y/z_mesh'
[values, values_3D] = gen_values(xsffile5, npts);
fprintf('Loaded in values from xsf')
fprintf('\n')
%Convert cartesian coordinates to spherical coordinates 'sphcoords'
[sphcoords] = gen_sphcoords(cartcoords);
fprintf('Converted from Cartesian to spherical coordinates')
fprintf('\n')
%Order the spherical coordinates to make it easier to integrate
[sphcoords_ordered,values_ordered]=order_sphscalarfield(sphcoords, values);
fprintf('Ordered spherical coordinates and corresponding values')
fprintf('\n')
%Limiting the 'index_cutoff'makes the code run way faster because the 
%interpolation step is super slow. For NiO, r<3 A at 99562 points
%Test: No real difference in F_k between r<3 and r<4 for NiO_00001 (dz^2)
index_cutoff = 99562;
sphcoords_ordered = sphcoords_ordered(1:index_cutoff,:);
values_ordered = values_ordered(1:index_cutoff,:);
%Integrate the wavefunction across theta and phi to get the radial wavefn R
%as a function of 'r_range'
fprintf('Beginning angular integration of wavefunction')
fprintf('\n')

[r_range5, R5] = integrate_angular(sphcoords_ordered,values_ordered);

%Take the average of all 3d orbital radial functions
R = (R1+R2+R3+R4+R5)/5;
%all r_range are the same
r = r_range1;
R0_dd = slater_integrate(r, R, r, R, 1, 0);
R2_dd = slater_integrate(r, R, r, R, 1, 2);
R4_dd = slater_integrate(r, R, r, R, 1, 4);



