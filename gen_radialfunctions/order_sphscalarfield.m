function [ sphcoords_ordered, values_ordered ] = order_sphscalarfield(...
    sphcoords, values)
%order_sphscalarfield : order 'sphcoords' in column-major mode and 'values'
%correspondingly to conduct multidimensional integration efficiently in 
%'slater_integrate'

%input: 'sphcoords' from 'gen_sphcoords', 'values' from 'gen_values'
%output: 'sphcoords_ordered' in column-major mode, 'values_ordered' ordered
%correspondingly

%Round to nearest hundredth because not rounding at all gives too much
%noise and differences between 'r' on the order of numerical error
%sphcoords = round(sphcoords,2);

%Merge 'sphcoords' and 'values' to sort them together.  
scalarfield = [sphcoords values];
%Use '[3 2 1]' input to sort the third column (r)in ascending order, and if
%there are ties, then sort by second (theta), then first (phi).
scalarfield_ordered = sortrows(scalarfield, [3 2 1]);
%split into 'sphcoords_ordered' and 'values_ordered'
sphcoords_ordered = scalarfield_ordered(:,[1 2 3]);
values_ordered = scalarfield_ordered(:,4);

end

