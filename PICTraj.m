 classdef PICTraj
  % PIC particle trajectory data
  %   pic = PIC(h5FilePath)
  
  properties
    id
    nt
    t
    x    
    y
    z
    vx
    vy
    vz
    Ex
    Ey
    Ez
    Bx
    By
    Bz
    t0
    x0
    y0
    z0
    vx0
    vy0
    vz0 
    mass
    charge
    userData = []; % anything can be added here
  end
  methods
    function obj = PICTraj(h5file)
      % traj = PICTraj('opt1',arg1,'opt2',arg2,...)
      %
      
%       if nargin > 0; have_input = 1; args = varargin; end
%       while have_input
%         l = 2;
%         switch lower(args{1})          
%           case {'file','h5file'}
%             h5file = args{2};
%             l = 2;
%           case ''
%         end
%         args = args{l+1:end};            
%       end
      
      if not(nargin == 0) % if its 0 then allocate the empty object array
        % get h5 info
        info = h5info(h5file);
        ntraj = numel(info.Groups(1).Groups); 
        iTr = read_itr(info);
        nTr = numel(iTr);
        obj(nTr,1) = obj;
        for itraj_ = 1:nTr
          itraj = iTr(itraj_);
          obj(itraj_).id = itraj;
          obj(itraj_).t = read_trajectories(itraj,1,'t');                              
          obj(itraj_).x = read_trajectories(itraj,1,'x');
          obj(itraj_).y = read_trajectories(itraj,1,'y');
          obj(itraj_).z = read_trajectories(itraj,1,'z');
          obj(itraj_).vx = read_trajectories(itraj,1,'vx');
          obj(itraj_).vy = read_trajectories(itraj,1,'vy');
          obj(itraj_).vz = read_trajectories(itraj,1,'vz');
          obj(itraj_).Ex = read_trajectories(itraj,1,'Ex');
          obj(itraj_).Ey = read_trajectories(itraj,1,'Ey');
          obj(itraj_).Ez = read_trajectories(itraj,1,'Ez');
          obj(itraj_).Bx = read_trajectories(itraj,1,'Bx');
          obj(itraj_).By = read_trajectories(itraj,1,'By');
          obj(itraj_).Bz = read_trajectories(itraj,1,'Bz');

          
          obj(itraj_).t0 = read_attributes(itraj,1,'t0');
          obj(itraj_).x0 = read_attributes(itraj,1,'x0');
          obj(itraj_).y0 = read_attributes(itraj,1,'y0');
          obj(itraj_).z0 = read_attributes(itraj,1,'z0');
          obj(itraj_).vx0 = read_attributes(itraj,1,'vx0');
          obj(itraj_).vy0 = read_attributes(itraj,1,'vy0');
          obj(itraj_).vz0 = read_attributes(itraj,1,'vz0');
          obj(itraj_).mass = read_attributes(itraj,1,'m');
          obj(itraj_).charge = read_attributes(itraj,1,'q');
          
          obj(itraj_) = obj(itraj_).rem_duplicates;
        end              
      end
      
      % these are unnecessary, can just read directly
      function out = read_itr(info)
        iGroup = find(contains({info.Groups.Name},'/traj'));
        nOutput = numel(info.Groups(iGroup).Groups);
        for iOutput = 1:nOutput
          str = info.Groups(iGroup).Groups(iOutput).Name;
          split_str = strsplit(str,'/');
          groups(iOutput) = str2num(split_str{3});
        end

        out = groups;      
      end
      function out = read_attributes(it,nt,att)
        % obj.read_attributes
        % get [x0,y0,z0,vx0,vy0,vz0,t0]

        init = zeros(nt,1);
        for it_ = 1:nt
          iTraj = it(it_);
          str_traj = sprintf('%06.0f',iTraj);
          init(it_,1) = h5readatt(h5file,['/traj/' str_traj],att);
        end     
        out = init;
      end 
      function out = read_trajectories(it,nt,field)
        % get indices      
        % initialize empty structure        
        for it_ = 1:nt
          iTraj = it(it_);
          str_traj = sprintf('%06.0f',iTraj);          
          out = h5read(h5file,['/traj/' str_traj '/' field]);%data_tmp; 
        end        
      end    
    end
    
    % Dervied quantities
    function out = U(obj)
      % PICTRAJ.U Kinetic energy of particle.
      
      for itr = 1:obj.ntr
        obj_tmp = obj(itr);
        out(itr).U = 0.5*obj_tmp.mass*(obj_tmp.vx.^2 + obj_tmp.vy.^2 + obj_tmp.vz.^2);
      end
      
      if obj.ntr == 1
        out = out.U;
      end
    end
    function out = Ux(obj)
      % PICTRAJ.U Kinetic energy of particle.
      
      for itr = 1:obj.ntr
        obj_tmp = obj(itr);
        out(itr).U = 0.5*obj_tmp.mass*(obj_tmp.vx.^2);
      end
      
      if obj.ntr == 1
        out = out.U;
      end
    end
    function out = Uy(obj)
      % PICTRAJ.U Kinetic energy of particle.
      
      for itr = 1:obj.ntr
        obj_tmp = obj(itr);
        out(itr).U = 0.5*obj_tmp.mass*(obj_tmp.vy.^2);
      end
      
      if obj.ntr == 1
        out = out.U;
      end
    end
    function out = Uz(obj)
      % PICTRAJ.U Kinetic energy of particle.
      
      for itr = 1:obj.ntr
        obj_tmp = obj(itr);
        out(itr).U = 0.5*obj_tmp.mass*(obj_tmp.vz.^2);
      end
      
      if obj.ntr == 1
        out = out.U;
      end
    end   
    function out = W(obj)
      % PICTRAJ.U Work done on the particle.
      % W = F dot dl = Fxdx + Fydy + Fzdz
      % F = qE
      
      for itr = 1:obj.ntr
        obj_tmp = obj(itr).rem_duplicates;
        dt = diff(obj_tmp.t);
        dx = diff(obj_tmp.x); %dx = interp1(obj_tmp.t(1:end-1)+dt,dx,obj_tmp.t); find(isnan(dx))
        dy = diff(obj_tmp.y); %dy = interp1(obj_tmp.t(1:end-1)+dt,dy,obj_tmp.t); find(isnan(dy))
        dz = diff(obj_tmp.z); %dz = interp1(obj_tmp.t(1:end-1)+dt,dz,obj_tmp.t); find(isnan(dz))
        out(itr).W = obj_tmp.charge*(obj_tmp.Ex.*[0;dx] + obj_tmp.Ey.*[0;dy] + obj_tmp.Ez.*[0;dz]);
      end
      
      if obj.ntr == 1
        out = out.W;
      end
    end
    function out = Wx(obj)
      % PICTRAJ.U Work done on the particle.
      % W = F dot dl = Fxdx + Fydy + Fzdz
      % F = qE
      
      for itr = 1:obj.ntr
        obj_tmp = obj(itr).rem_duplicates;
        dt = diff(obj_tmp.t);
        dx = diff(obj_tmp.x); %dx = interp1(obj_tmp.t(1:end-1)+dt,dx,obj_tmp.t); find(isnan(dx))
        
        out(itr).W = obj_tmp.charge*(obj_tmp.Ex.*[0;dx]);
      end
      
      if obj.ntr == 1
        out = out.W;
      end
    end
    function out = Wy(obj)
      % PICTRAJ.U Work done on the particle.
      % W = F dot dl = Fxdx + Fydy + Fzdz
      % F = qE
      
      for itr = 1:obj.ntr
        obj_tmp = obj(itr).rem_duplicates;
        dt = diff(obj_tmp.t);
        dx = diff(obj_tmp.x); %dx = interp1(obj_tmp.t(1:end-1)+dt,dx,obj_tmp.t); find(isnan(dx))
        dy = diff(obj_tmp.y); %dy = interp1(obj_tmp.t(1:end-1)+dt,dy,obj_tmp.t); find(isnan(dy))
        dz = diff(obj_tmp.z); %dz = interp1(obj_tmp.t(1:end-1)+dt,dz,obj_tmp.t); find(isnan(dz))
        out(itr).W = obj_tmp.charge*(obj_tmp.Ey.*[0;dy]);
      end
      
      if obj.ntr == 1
        out = out.W;
      end
    end
    function out = Wz(obj)
      % PICTRAJ.U Work done on the particle.
      % W = F dot dl = Fxdx + Fydy + Fzdz
      % F = qE
      
      for itr = 1:obj.ntr
        obj_tmp = obj(itr).rem_duplicates;
        dt = diff(obj_tmp.t);
        dx = diff(obj_tmp.x); %dx = interp1(obj_tmp.t(1:end-1)+dt,dx,obj_tmp.t); find(isnan(dx))
        dy = diff(obj_tmp.y); %dy = interp1(obj_tmp.t(1:end-1)+dt,dy,obj_tmp.t); find(isnan(dy))
        dz = diff(obj_tmp.z); %dz = interp1(obj_tmp.t(1:end-1)+dt,dz,obj_tmp.t); find(isnan(dz))
        out(itr).W = obj_tmp.charge*(obj_tmp.Ez.*[0;dz]);
      end
      
      if obj.ntr == 1
        out = out.W;
      end
    end 
    % Analysis
    function out = zcross(obj)
      % PICTraj.ZCROSS Locations where the particle crossez z = 0.
      % See also PICTRAJ.NCROSS
      for itr = 1:obj.ntr
        t = obj(itr).t;
        z = obj(itr).z;
        x = obj(itr).x;
        %idup = find(diff(obj(itr).t)==0);
        %idup = find(diff(obj(itr).z)<1e-11);
        %t(idup) = [];
        %x(idup) = [];
        %z(idup) = [];
                
        zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
        icross = zci(z);
        icross(icross == 1) = [];
        icross(icross == numel(z)) = [];
        
        idup = find(diff(icross)==1);
        icross(idup) = [];
        
        out(itr).ic = icross;
        out(itr).t = t(icross);
        out(itr).x = x(icross);
        out(itr).z = z(icross);
        out(itr).nc = numel(icross);
                
        if 1 % plot
          %plot(obj(itr).x,obj(itr).z,obj(itr).x(icross),obj(itr).z(icross),'*')
          plot(x,z,x(icross),z(icross),'*')
          set(gca,'XGrid','on','YGrid','on')
          title(gca,sprintf('itr = %.0f, ncr = %g',itr,out(itr).nc))
          pause
        end
      end      
    end
    function out = ncross(obj)
      % PICTraj.NCROSS Number of crosses where the particle crossez z = 0.
      % See also PICTRAJ.ZCROSS
      for itr = 1:obj.ntr
        z = obj(itr).z;
        zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
        icross = zci(z);
        icross(icross == 1) = [];
        icross(icross == numel(z)) = [];
        idup = find(diff(icross)==1);
        icross(idup) = [];
        
        xcr_all(itr) = numel(icross);
        if 0 % plot
          plot(obj(itr).x,obj(itr).z,obj(itr).x(icross),obj(itr).z(icross),'*')
          set(gca,'XGrid','on','YGrid','on')
          title(hca,sprintf('itr = %.0f',itr))
          pause
        end
      end
      out = xcr_all;
    end
    
    % Operators
    function out = nt_(obj)
      value = numel(obj);
      out = value;
    end
    function out = length(obj)
      for itr = 1:obj.ntr
        value(itr) = numel(obj(itr).t);
      end
      out = value;
    end
    function out = max(obj,field)
      maxval = nan(numel(obj),1);
      for itr = 1:numel(obj)
        
      end
      
    end
    function out = nancat(obj)
      % PICTRAJ.NANCAT Don't know what it was supposed to do, not
      % implemented.
      %TR(itr)
      out = obj;
    end
    %
    function TR = rem_duplicates(obj)
    % PICTRAJ.REM_DUPLICATES Removes duplicates points, which may cause
    %   problems when interpolating data;
    % 
      for itr = 1:obj.ntr
        idup = find(diff(obj(itr).t)==0);
        ikeep = find(not(diff(obj(itr).t)==0));
        ikeep = setdiff(1:obj(itr).length,idup);
        obj(itr) = obj(itr).select_inds(ikeep);
      end
      TR = obj;
    end
    % Finding subsets of data
    function TR = find(obj,varargin)
      TR = obj;%(find(obj));
    end
    function TR = pass(obj,varargin)
      if mod(nargin,1) == 1
        error('Wrong number of input.')
      end
      
      ikeep = [];
      
      for itr = 1:obj.ntr
        ind = 1:numel(obj(itr).t); % at beginning all inds are ok
        args = varargin;    
        nargs = numel(args);    
        for iarg = 1:(nargs/2) % switch for while, incase conditions is met early.
          field = args{1};
          lim = args{2};
          args = args(3:end);
          
          if isfield(obj,field)    % is field
            data = obj(itr).(field);            
          else % is expression
            % excepted expressions
            exceptions = {'atan','cosd','angle','sind'};
            replacements = {'$','#','@','&'};
            field = replace(field,exceptions,replacements);
            newfields = cellfun(@(x) sprintf('obj(itr).%s',x),obj.fieldnames,'UniformOutput',0);
            newexp = replace(field,obj.fieldnames,newfields);
            newexp = replace(newexp,replacements,exceptions);
            data = eval(newexp); 
          end
          ind = intersect(ind,find(data>=lim(1)));
          ind = intersect(ind,find(data<=lim(2)));
        end
        if not(isempty(ind))
          %disp(sprintf('keeping %g',itr))
          ikeep = [ikeep itr];
        end
      end
      TR = obj(ikeep);      
      function out = replace_field_with_expression(fields)
        for ifield = 1:numel(fields)
          strcell{ifield} = sprintf('obj(itr).(%s)',fields{ifield}); 
        end
      end
    end
    function TR = lim(obj,varargin)
      if mod(nargin,1) == 1
        error('Wrong number of input.')
      end
      
      ikeep = [];
      
      TR = obj;
      
      ikeep = [];
      for itr = 1:obj.ntr
        %disp(sprintf('iTr = %g',itr))
        ind = 1:numel(obj(itr).t); % at beginning all inds are ok
        args = varargin;    
        nargs = numel(args);  
        
        for iarg = 1:(nargs/2) % switch for while, incase conditions is met early.
          field = args{1};
          lim = args{2};
          args = args(3:end);
          
          data = obj(itr).(field);
          ind = intersect(ind,find(data>=lim(1)));
          ind = intersect(ind,find(data<=lim(2)));          
        end
        if not(isempty(ind))
          ikeep = [ikeep itr];
        end
        TR(itr) = TR(itr).select_inds(ind);  
      end   
      TR = TR(ikeep);   
    end    
    function TR = select_inds(obj,inds)
      % PICTraj
      TR = obj;
      ndata = numel(TR.t);
      fields = fieldnames(obj);
      
      for ifield = 1:numel(fields)
        data = TR.(fields{ifield});
        if not(numel(data)==ndata)
          continue
        end
        TR.(fields{ifield}) = data(inds);
      end
    end
    
    % 
    function value = ntr(obj)
      value = numel(obj);
    end
    
    % Plotting
    function varargout = plot_single(obj,varstrs,varargin)
      % PICTRAJ.PLOT_SINGLE Plots single trajectory in various ways.
      %   PICTRAJ.PLOT_SINGLE('opt1','opt2',...'optN')
      %   Possible options:
      %   any combination of x y z - trajectory in given plane, or 3D
      %   any combination of vx vy vz - velocity in phase space
            % Default options, values
      doAdjustCLim = 0;
      cmap = pic_colors('blue_red');
      doAdjustCMap = 0;
      nvars = numel(varstrs);
      
      have_options = 0;
      nargs = numel(varargin);      
      if nargs > 0, have_options = 1; args = varargin(:); end
      
      while have_options
        l = 1;
        varstr = lower(args{1});
        switch varstr
          case {'xy','yz','xz','xz','yz','zy'}
            do2Dtraj = 1;
            coord2D = args{2};
            l = 2;            
          case {'xyz','yzx','zxy','zyx','yxz','xzy'}
            do3Dtraj = 1;
            coord3D = args{2};
            l = 2;
          case {'vxvy','vyvz','vxvz','vxvz','vyvz','vzvy'}
            doV2Dtraj = 1;
            coordV2D = args{2};
            l = 2;
          case 'clim'
            l = 2;
            doAdjustCLim = 1;  
            clims = args{2};
          case 'cmap'
            l = 2;
            doAdjustCMap = 1;
            cmaps = args{2};        
          otherwise
            warning(sprintf('Input ''%s'' not recognized.',args{1}))            
        end        
        args = args(l+1:end);
        if isempty(args), break, end    
      end
                  
      % setup figure
      fig = figure;      
      [nrows,ncols] = size(varstrs);           
      npanels = nrows*ncols;
      ip = 0;
      for irow = 1:nrows
        for icol = 1:ncols
          ip = ip + 1;
          h(irow,icol) = subplot(nrows,ncols,ip);
        end
      end
      
      for ivar = 1:nvars
        hca = h(ivar);
        varstr_split = strsplit(varstrs{ivar},'_');
        varstr = varstr_split{1};
        switch varstr
          case {'tU','tW'}
            plot(hca,obj.(varstr(1)),obj.(varstr(2)))
            hca.XLabel.String = sprintf('%s ()',varstr(1));
            hca.YLabel.String = sprintf('%s ()',varstr(2)); 
            hca.XGrid = 'on';
            hca.YGrid = 'on';
          case {'xU','yU','zU','xUx','yUx','zUx','xUy','yUy','zUy','xUz','yUz','zUz',...
              'xW','yW','zW','xWx','yWx','zWx','xWy','yWy','zWy','xWz','yWz','zWz',...
              'xEx','yEx','zEx','xEy','yEy','zEy','xEz','yEz','zEz',...
              'xvx','yvx','zvx','xvy','yvy','zvy','xvz','yvz','zvz'}
            if numel(varstr_split) > 1 && strcmp(varstr_split{2},'cumsum')
              plot(hca,obj.(varstr(1)),cumsum(obj.(varstr(2:end)),'omitnan'))
            else
              plot(hca,obj.(varstr(1)),obj.(varstr(2:end)))
            end            
            hca.XLabel.String = sprintf('%s ()',varstr(1));
            hca.YLabel.String = sprintf('%s ()',varstr(2:end)); 
            hca.XGrid = 'on';
            hca.YGrid = 'on';
          case {'xy','yz','xz','xz','yz','zy'}
            plot(hca,obj.(varstr(1)),obj.(varstr(2)))
            hca.XLabel.String = sprintf('%s (d_i)',varstr(1));
            hca.YLabel.String = sprintf('%s (d_i)',varstr(2)); 
            hca.XGrid = 'on';
            hca.YGrid = 'on';
          case {'xyz','yzx','zxy','zyx','yxz','xzy'}
            plot2(hca,obj.(varstr(1)),obj.(varstr(2)),obj.(varstr(3)))
            hca.XLabel.String = sprintf('%s (d_i)',varstr(1));
            hca.YLabel.String = sprintf('%s (d_i)',varstr(2)); 
            hca.YLabel.String = sprintf('%s (d_i)',varstr(3)); 
            hca.XGrid = 'on';
            hca.YGrid = 'on';
          case {'vxvy','vyvz','vxvz','vxvz','vyvz','vzvy',...
                'vxEx','vxEy','vxEz','vyEx','vyEy','vyEz','vzEx','vzEy','vzEz'}
            plot(hca,obj.(varstr(1:2)),obj.(varstr(3:4)));
            hca.XLabel.String = sprintf('%s (d_i)',varstr(1:2));
            hca.YLabel.String = sprintf('%s (d_i)',varstr(3:4)); 
            hca.XGrid = 'on';
            hca.YGrid = 'on';
          otherwise
            warning(sprintf('Variable %s not supported/implemented.',varstr_split{1}))
        end
      end
    end
    function h = plot_all_xz(obj,varargin)
      % PICTRAJ.PLOT_ALL_XZ Plots all trajectories in xz plane.
      
      for itr = 1:obj.ntr
        plot(obj(itr).x,obj(itr).z);
        if itr == 1
          hold(gca,'on')
        end
      end
      hold(gca,'off')
      h=gca;      
      h.XGrid = 'on';
      h.YGrid = 'on';
      h.XLabel.String = 'x (d_i)';
      h.YLabel.String = 'z (d_i)';
    end
    function h = plot_all_xy(obj,varargin)
      % PICTRAJ.PLOT_ALL_XY Plots all trajectories in xy plane.
      
      for itr = 1:obj.ntr
        plot(obj(itr).x,obj(itr).y);
        if itr == 1
          hold(gca,'on')
        end
      end
      hold(gca,'off')
      h=gca;      
      h.XGrid = 'on';
      h.YGrid = 'on';
      h.XLabel.String = 'x (d_i)';
      h.YLabel.String = 'y (d_i)';
    end
    function h = plot_all_xyz(obj,varargin)
      % PICTRAJ.PLOT_ALL_XZZ Plots all trajectories in xyz plane.
      
      for itr = 1:obj.ntr
        plot3(obj(itr).x,obj(itr).y,obj(itr).z);
        if itr == 1;
          hold(gca,'on')
        end
      end
      hold(gca,'off')
      h=gca;
      view(h,[0 1 0])
      h.XGrid = 'on';
      h.YGrid = 'on';
      h.ZGrid = 'on';
      h.XLabel.String = 'x (d_i)';
      h.YLabel.String = 'y (d_i)';
      h.ZLabel.String = 'z (d_i)';
    end
    
    %
    function out = fieldnames(obj)
      out = properties(obj);
    end
    function out = isfield(obj,field)
      if find(strcmp(obj.fieldnames,field))
        out = 1;
      else
        out = 0;
      end        
    end
  end
  
end