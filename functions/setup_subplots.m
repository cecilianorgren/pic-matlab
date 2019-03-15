function h = setup_subplots(nrows,ncols,orientation)

orientation_default = 'horizontal';
if not(exist('orientation','var'))
  orientation = orientation_default;
elseif not(any(strcmp(orientation,{'horizontal','vertical'})))
  warning(sprintf('Orientation: %s, not recognized. Using default orientation %s',orientation,orientation_default))
  orientation = orientation_default;
end
  
npanels = nrows*ncols;
ipanel = 0;
if strcmp(orientation,'vertical')
  for irow = 1:nrows
    for icol = 1:ncols
      ipanel = ipanel + 1;
      h(irow,icol) = subplot(nrows,ncols,ipanel);
      h_ind(irow,icol) = ipanel;
    end
  end  
else
  for icol = 1:ncols
    for irow = 1:nrows
      ipanel = ipanel + 1;
      h(irow,icol) = subplot(nrows,ncols,ipanel);
      h_ind(irow,icol) = ipanel;
    end
  end
end
%%
%if strcmp(orientation,'vertical')
%  h = h';
%end
h = h(:);