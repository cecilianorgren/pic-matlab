function out = track_xlines(input)

dx_max = 10; % max change in dx (grid cells) for it to be considered the same xline
sim = input;
x = sim.xi;
z = sim.zi;
xlines = cell(0,0);
for it = 1:sim.length
  disp(sprintf('%g/%g',it,sim.length));
  % load B
  Bx = squeeze(sim(it).Bx);
  Bz = squeeze(sim(it).Bz);
  A = vector_potential(x,z,Bx,Bz); % vector potential
  [saddle_locations,saddle_values] = saddle(A,'sort');  
  if it == 1
    for isaddle = 1:numel(saddle_values)
      xlines{isaddle} = {it,loc,val};            
    end
  else
    nold = numel(xlines);
    nnew = numel(saddle_values);
    xoldnew = nchoosek(1:nold,nnew); % check all these pairs to see if they belong together
  end
  if not(isempty(saddle_locations)) 
    for isaddle = 1:numel(saddle_values)
      if it == 1
        xlines{isaddle} = {it,loc,val};
        continue
      end      
      % can only match one existing xline with one new saddle location
      
      
    end
  end
%   if not(isempty(saddle_locations)) 
%     for isaddle = 1:numel(saddle_values)
%       % can only match one existing xline with one new saddle location
%       xlines = check_if_new_ord_existing(xlines,saddle_locations(isaddle,:),saddle_values(isaddle));
%       for iix = 1:numel(xlines); plot(xlines{iix}{1},xlines{iix}{2}(:,1),'.-'); drawnow; hold on; end
%     end
%   end
end
out = xlines;

function xl = check_if_new_ord_existing(xl,loc,val)

  found_xline = 0;
  if isempty(xl) % first iteration
    xl{1} = {it,loc,val};
    found_xline = 1;    
  else
    % run through existing xlines
    for ix = 1:numel(xl)
      dx = loc(1) -  xl{ix}{2}(1);
      dx
      if abs(dx) < dx_max % this xline corresponds to existing xline
        old_t = xl{ix}{1};
        old_xy = xl{ix}{2};
        old_val = xl{ix}{3};
        xl{ix} = {[old_t;it],[old_xy;loc],[old_val;val]};
        found_xline = 1;
      end
    end
    % didnt find existing xline, so add new
    if not(found_xline)
      xl{end+1} = {it,loc,val};
    end
  end
end

end