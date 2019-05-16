function out = interpolate_time(old_points,old_data,new_points)


n_new = numel(new_points);
[n_old,nx,nz] = size(old_data);

new_data = nan(n_new,nx,nz);


% display_progress at
ix_display = round(nx/nx):round(nx/nx):nx;
ix_display = 100:100:nx;
for ix = 1:nx
  if mod(ix,ix_display(1)) == 0
    ix_display = ix_display(2:end);
    disp(sprintf('ix = %g/%g',ix,nx))
  end
  for iz = 1:nz
    if 1 % interp1 function
      new_data(:,ix,iz) = interp1(old_points,old_data(:,ix,iz),new_points);
    else % simple linear interpolation
      %dt_plus = ;
      new_data(:,ix,iz) = interp1(old_points,old_data(:,ix,iz),new_points);
    end
  end
end

out = new_data;