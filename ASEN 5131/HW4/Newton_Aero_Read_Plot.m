% -- Reduced Newton_Aero.m to demonstrate reading and plotting

clc; clear all; close all;

grid_fname = 'Sphere_Processed.dat';

% ===== Read and Process Grid
% -- One row per node, first 3 nodes form 1st triangle, and so on
grid_raw = dlmread(grid_fname);

n_tri = size(grid_raw,1)/3;
if (mod(n_tri,1)~=0)
   error('Error, number of nodes is not divisible by 3\n');
end

% -- grid_x(k,j,i) has k-th dimension of j-th node of surf tri i
grid_x = zeros(3,3,n_tri);
center_x = zeros(3,n_tri);
for i=1:n_tri
   for j=1:3
      grid_x(1:3,j,i) = grid_raw( 3*(i-1)+j, 1:3);
      center_x(:,i) = center_x(:,i) + grid_x(:,j,i)/3;
   end
end

% ===== Plot coefficient of pressure on triangulated surface
% -- Fill C_p with dummy data, pick z of triangle center plus 1
C_p = center_x(3,:)+1;

% -- Node indexes that define each triangle: (1,2,3), (4,5,6), ...
% -- If grid format was more sophisticated than STL, this would change
tri_conn = zeros(n_tri,3);
for i=1:n_tri
   tri_conn(i,1:3) = [ 3*(i-1)+1, 3*(i-1)+2, 3*(i-1)+3 ];
end

% -- Plot and label
trisurf(tri_conn, grid_raw(:,1), grid_raw(:,2), grid_raw(:,3), C_p);
title('Cp');
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
colorbar; caxis([ 0 2 ]);
