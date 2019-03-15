function handles = draw_circle(x,z,circle_radius,circle_axis)
% DRAW_CIRCLE draws circle with radius circle_radius and axis circle_axis.
%   handles = draw_circle(x,z,circle_radius,circle_axis)

nangles = 30;
angle = linspace(0,2*pi,nangles);

ncircles = numel(circle_radius);

r1 = circle_axis; % r2 and r3 are interchangeable
r2 = cross_product(r1.x,r1.y,r1.z,0,1,0);
r2.abs = sqrt(r2.x.^2 + r2.y.^2 + r2.z.^2);
r2.x = r2.x./r2.abs;
r2.y = r2.y./r2.abs;
r2.z = r2.z./r2.abs;
r2.abs = sqrt(r2.x.^2 + r2.y.^2 + r2.z.^2);
r2 = cross_product(r2.x,r2.y,r2.z,r1.x,r1.y,r1.z);
r3 = cross_product(r1.x,r1.y,r1.z,r2.x,r2.y,r2.z);
r3.abs = sqrt(r3.x.^2 + r3.y.^2 + r3.z.^2);
  
ax = gca;
hold(ax,'on')
for icircle = 1:ncircles
  % circle with axis along x
  xx = 0;
  yy = circle_radius(icircle)*cos(angle);
  zz = circle_radius(icircle)*sin(angle);
  
  if not(all(isnan(xx))) && x(icircle)<70 && x(icircle)>0 && z(icircle)<8 && z(icircle)>-8
    r1_tmp = [r1.x(icircle) r1.y(icircle) r1.z(icircle)];
    r2_tmp = [r2.x(icircle) r2.y(icircle) r2.z(icircle)];
    r3_tmp = [r3.x(icircle) r3.y(icircle) r3.z(icircle)];
    

    % Rotate 
    [new_x,new_y,new_z] = rotate_xyz(xx,yy,zz,r1_tmp,r2_tmp,r3_tmp);  
    plot(ax,x(icircle)+new_x,z(icircle)+new_z)
    drawnow;
  end
end
hold(ax,'off')  