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
    A
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
      %tic; % time it to get a feel for how different file sizes slow affect the loading
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
          obj(itraj_).A = read_trajectories(itraj,1,'Ay');

          
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
      %toc
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
    function out = Ugen(obj,comp,ind)
      % PICTRAJ.U Kinetic energy of particle.
      % out = Ugen(obj,comp,ind)
      %
        
      if or(not(exist('ind','var')),isempty(ind))       
        istart = repmat(1,obj.ntr,1);
        istop = [obj.length];
      elseif numel(ind) == 1
        istart = repmat(ind,obj.ntr,1);
        istop = istart;
      elseif numel(ind) == obj.ntr
        istart = ind;
        istop = istart;
      else
        error(sprintf('Wrong dimension of input.'))
      end
      
      switch comp
        case {'x','y','z'}
          for itr = 1:obj.ntr
            obj_tmp = obj(itr);
            U_tmp = 0.5*obj_tmp.mass*obj_tmp.(['v' comp]).^2;
            out(itr).U = U_tmp(istart(itr):istop(itr));
          end
        case {'xyz','tot'}
          for itr = 1:obj.ntr
            obj_tmp = obj(itr);  
            U_tmp = 0.5*obj_tmp.mass*(obj_tmp.vx.^2 + obj_tmp.vy.^2 + obj_tmp.vz.^2);
            out(itr).U = U_tmp(istart(itr):istop(itr));
          end        
      end      
      if obj.ntr == 1
        out = out.U;
      elseif istart == istop
        out = [out.U];
      end
    end
    function out = U(obj)
      % PICTRAJ.U Kinetic energy of particle.      
      out = Ugen(obj,'tot',[]);
    end
    function out = Ux(obj)
      % PICTRAJ.Ux Kinetic energy of particle.
      out = Ugen(obj,'x',[]);
    end
    function out = Uy(obj)
      % PICTRAJ.Uy Kinetic energy of particle.
      out = Ugen(obj,'y',[]);
    end
    function out = Uz(obj)
      % PICTRAJ.Uz Kinetic energy of particle.
      out = Ugen(obj,'z',[]);      
    end   
    function out = Ustart(obj)
      % PICTRAJ.U Kinetic energy of particle.
      out = Ugen(obj,'tot',1);    
    end  
    function out = Uxstart(obj)
      % PICTRAJ.U Kinetic energy of particle.
      out = Ugen(obj,'x',1);    
    end  
    function out = Uystart(obj)
      % PICTRAJ.U Kinetic energy of particle.      
      out = Ugen(obj,'y',1);
    end  
    function out = Uzstart(obj)
      % PICTRAJ.U Kinetic energy of particle.      
      out = Ugen(obj,'z',1);     
    end  
    function out = Ustop(obj)
      % PICTRAJ.U Kinetic energy of particle.      
      out = Ugen(obj,'tot',obj.length);     
    end  
    function out = Uxstop(obj)
      % PICTRAJ.U Kinetic energy of particle.      
      out = Ugen(obj,'x',obj.length);
    end  
    function out = Uystop(obj)
      % PICTRAJ.U Kinetic energy of particle.      
      out = Ugen(obj,'y',obj.length);
    end  
    function out = Uzstop(obj)
      % PICTRAJ.U Kinetic energy of particle.      
      out = Ugen(obj,'z',obj.length);
    end
    function out = dU(obj)
      out = obj.Ustop-obj.Ustart;
    end
    function out = dUz(obj)
      out = obj.Uzstop-obj.Uzstart;
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
        dy = diff(obj_tmp.y); %dy = interp1(obj_tmp.t(1:end-1)+dt,dy,obj_tmp.t); find(isnan(dy))
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
        dz = diff(obj_tmp.z); %dz = interp1(obj_tmp.t(1:end-1)+dt,dz,obj_tmp.t); find(isnan(dz))
        out(itr).W = obj_tmp.charge*(obj_tmp.Ez.*[0;dz]);
      end
      
      if obj.ntr == 1
        out = out.W;
      end
    end 
    function out = Wsum(obj)
      for itr = 1:obj.ntr
        obj_tmp = obj(itr).rem_duplicates;        
        W = obj_tmp.W;
        out(itr) = sum(W);
      end
    end
    function out = Wxsum(obj)
      for itr = 1:obj.ntr
        obj_tmp = obj(itr).rem_duplicates;        
        W = obj_tmp.Wx;
        out(itr) = sum(W);
      end
    end
    function out = Wysum(obj)
      for itr = 1:obj.ntr
        obj_tmp = obj(itr).rem_duplicates;        
        W = obj_tmp.Wy;
        out(itr) = sum(W);
      end
    end
    function out = Wzsum(obj)
      for itr = 1:obj.ntr
        obj_tmp = obj(itr).rem_duplicates;        
        W = obj_tmp.Wz;
        out(itr) = sum(W);
      end
    end
    function out = vB(obj,comp)
      % any combination of vB
      % out = vB(obj,'xy')
      for itr = 1:obj.ntr
        obj_tmp = obj(itr).rem_duplicates;                  
        v = obj_tmp.(['v' comp(1)]);
        B = obj_tmp.(['B' comp(2)]);                  
        out(itr).vB = v.*B;
      end
      
      if obj.ntr == 1
        out = out.vB;
      end
    end
    function out = vxBy(obj)
      % See also PICTraj.vB
      out = obj.vB('xy');
    end
    function out = vyBx(obj)
      % See also PICTraj.vB
      out = obj.vB('yx');
    end
    function out = vxBz(obj)
      % See also PICTraj.vB
      out = obj.vB('xz');
    end
    function out = vzBx(obj)
      % See also PICTraj.vB
      out = obj.vB('zx');
    end
    function out = vyBz(obj)
      % See also PICTraj.vB
      out = obj.vB('yz');
    end
    function out = vzBy(obj)
      % See also PICTraj.vB
      out = obj.vB('zy');
    end
    function out = py(obj)
      %Ay = obj.A;
      
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
                
        if 0 % plot
          %plot(obj(itr).x,obj(itr).z,obj(itr).x(icross),obj(itr).z(icross),'*')
          plot(x,z,x(icross),z(icross),'*')
          set(gca,'XGrid','on','YGrid','on')
          title(gca,sprintf('itr = %.0f, ncr = %g',itr,out(itr).nc))
          pause
        end
      end      
    end
    function out = firstcross(obj)
      allcross = obj.zcross;
      for itr = 1:obj.ntr
        out(itr) = allcross(itr).x(1);
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
    
    % 
    function out = coordgen(obj,comp,ind)
      % PICTRAJ.XGEN x position energy of particle.
      % out = Ugen(obj,comp,ind)
      %
        
      if or(not(exist('ind','var')),isempty(ind))       
        istart = repmat(1,obj.ntr,1);
        istop = [obj.length];
      elseif numel(ind) == 1
        istart = repmat(ind,obj.ntr,1);
        istop = istart;
      elseif numel(ind) == obj.ntr
        istart = ind;
        istop = istart;
      else
        error(sprintf('Wrong dimension of input.'))
      end
      
      for itr = 1:obj.ntr
        obj_tmp = obj(itr);
        tmp = obj_tmp.(comp);
        out(itr).x = tmp(istart(itr):istop(itr));
      end                
      if obj.ntr == 1
        out = out.x;
      elseif istart == istop
        out = [out.x];
      end
    end
    function out = tstart(obj)
      out = obj.coordgen('t',1);      
    end
    function out = tstop(obj)
      out = obj.coordgen('t',obj.length);      
    end
    function out = xstart(obj)
      out = obj.coordgen('x',1);      
    end
    function out = xstop(obj)
      out = obj.coordgen('x',obj.length);      
    end
    function out = ystart(obj)
      out = obj.coordgen('y',1);      
    end
    function out = ystop(obj)
      out = obj.coordgen('y',obj.length);      
    end
    function out = zstart(obj)
      out = obj.coordgen('z',1);      
    end
    function out = zstop(obj)
      out = obj.coordgen('z',obj.length);      
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
      % PICTARAJ.FIND Finds subsets of data
      %   PICTARAJ.FIND(cond1,cond1,cond3,cond4)
      %   PICTARAJ.FIND([tr.t0]==160,[tr.x0]<170,tr.ncross==1)
      inds = 1:obj.ntr;
      for iarg = 1:numel(varargin)
        inds = intersect(inds,find(varargin{iarg}));
      end      
      TR = obj(unique(inds));
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
          
          if isfield(obj,field) % is field
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
    function TR = tlim(obj,tint)
      TR = obj.lim('t',tint); 
      for itr = 1:TR.ntr
        TR(itr) = TR(itr).lim('t',tint);
      end
    end
    function TR = subset(obj,inds)
      % PICTRAJ.SUBSET Same as PICTRAJ.SELECT_INDS
      % See also PICTRAJ.SELECT_INDS
      TR = obj.select_inds(inds);
    end
    function TR = select_inds(obj,inds)
      % PICTRAJ.SELECT_INDS Select subset of individual trajectory.
      % TR = TR.SELECT_SUBSET(1:10)
      TR = obj;      
      fields = fieldnames(TR);
      if numel(inds) == 1
        inds = repmat(inds,TR.ntr,1);
      elseif numel(inds) == TR.ntr
        inds = reshape(inds,numel(inds),1);
      elseif any(size(inds)==1)
        inds = repmat(reshape(inds,1,numel(inds)),TR.ntr,1);
      end
      
      for itr = 1:TR.ntr
        nt = TR(itr).length;
        for ifield = 1:numel(fields)
          data = TR(itr).(fields{ifield});
          if not(numel(data) == nt)
            continue
          end
          TR(itr).(fields{ifield}) = data(inds(itr,:));
        end
      end
    end
    function out = interp(obj,field,time)
      % PICTRAJ.INTERP Interpolates trjectories to given time/timeline.
      
      for itr = 1:obj.ntr
        tt = obj(itr).t;
        dd = obj(itr).(field);
        idup = find(diff(tt)==0);
        tt(idup) = [];
        dd(idup) = [];        
        out(itr).data = interp1(tt,dd,time,'linear','extrap');        
      end
      if numel(time) == 1
        out = [out.data];
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
      %fig = figure;      
      [nrows,ncols] = size(varstrs);           
      npanels = nrows*ncols;
      ip = 0;
      for irow = 1:nrows
        for icol = 1:ncols
          ip = ip + 1;
          h(irow,icol) = subplot(nrows,ncols,ip);
        end
      end
      
      %for itr = 1:obj.ntr
      for ivar = 1:nvars
        hca = h(ivar);
        varstr_split_vars = strsplit(varstrs{ivar},',');
        nvars = numel(varstr_split_vars);
        legs = cell(0);
        for ivar = 1:nvars
          varstr_split = strsplit(varstr_split_vars{ivar},'_');
          varstr = varstr_split{1};
          holdOn = 0;
          switch varstr
            case {'tU','tUx','tUy','tUz',...
                'tW','tWx','tWy','tWz',...
                'tEx','tEy','tEz',...
                'tBx','tBy','tBz',...
                'tEx','tEy','tEz',...
                'tvxBy','tvyBx','tvxBz','tvzBx','tvyBz','tvzBy',...
                'tvx','tvy','tvz',...
                'tx','ty','tz'}
              if numel(varstr_split) > 1 && strcmp(varstr_split{2},'cumsum')
                plot(hca,obj.(varstr(1)),cumsum(obj.(varstr(2:end)),'omitnan'))
              else
                plot(hca,obj.(varstr(1)),obj.(varstr(2:end)))
              end                        
              %plot(hca,obj.(varstr(1)),obj.(varstr(2:end)))
              hca.XLabel.String = sprintf('%s ()',varstr(1));
              hca.YLabel.String = sprintf('%s ()',varstr(2:end)); 
              hca.XGrid = 'on';
              hca.YGrid = 'on';
            case {'xU','yU','zU','xUx','yUx','zUx','xUy','yUy','zUy','xUz','yUz','zUz',...
                'xW','yW','zW','xWx','yWx','zWx','xWy','yWy','zWy','xWz','yWz','zWz',...
                'xEx','yEx','zEx','xEy','yEy','zEy','xEz','yEz','zEz',...
                'xBx','yBx','zBx','xBy','yBy','zBy','xBz','yBz','zBz',...
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
            case {'xy','yz','xz','yx','zy','zx'}
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
            case {'vxvy','vxvz','vyvx','vzvx','vyvz','vzvy',...
                  'vxEx','vxEy','vxEz','vyEx','vyEy','vyEz','vzEx','vzEy','vzEz',...
                  'WxWy','WxWz','WyWx','WzWx','WyWz','WzWy'}
              plot(hca,obj.(varstr(1:2)),obj.(varstr(3:4)));
              hca.XLabel.String = sprintf('%s (d_i)',varstr(1:2));
              hca.YLabel.String = sprintf('%s (d_i)',varstr(3:4)); 
              hca.XGrid = 'on';
              hca.YGrid = 'on';
            otherwise
              warning(sprintf('Variable %s not supported/implemented.',varstr_split{1}))          
          end
          legs{end+1} = varstr_split_vars{ivar};
          if ivar == 1 && holdOn  == 0
            hold(hca,'on')
            holdOn = 1;
          end
        end
        
        if nvars > 1
          legend(hca,legs,'location','west','interpreter','none')
        end
      end
      if obj.ntr == 1
        h(1).Title.String = sprintf('t_0w_{ci} = %.2f, [x_0,y_0,z_0] = [%.2f,%.2f,%.2f], [v_{x0},v_{y0},v_{z0}] = [%.2f,%.2f,%.2f]',obj.t0,obj.x0,obj.y0,obj.z0,obj.vx0,obj.vy0,obj.vz0);
      end
    end
    function h = plot_all(obj,xstr,ystr,varargin)
      % PICTRAJ.PLOT_ALL_XY Plots all trajectories in xy plane.
      
      for itr = 1:obj.ntr
        plot(obj(itr).(xstr),obj(itr).(ystr));
        if itr == 1
          hold(gca,'on')
        end
      end
      hold(gca,'off')
      h=gca;      
      h.XGrid = 'on';
      h.YGrid = 'on';
      h.XLabel.String = sprintf('%s (d_i)',xstr);
      h.YLabel.String = sprintf('%s (d_i)',ystr);
    end
    function h = plot_all_xz(obj,varargin)
      % PICTRAJ.PLOT_ALL_XZ Plots all trajectories in xz plane.
      
      % Defaults
      color = [0 0 0];
      doColor = 1;
      cmap = pic_colors('waterfall'); % maybe I should build-in the cmap
      have_input = 0;
      
      [h,args,nargs] = axescheck(varargin{:});
      if isempty(h)
        h = gca; % current or new axes
      end
      
      if nargs > 0; have_input = 1; args = args; end
      while have_input
        l = 2;
        switch lower(args{1})
          case {'color'}
            doColor = 1;
            color = args{2};
            l = 2;
          otherwise
            warning(sprintf('Unknown input %s.',args{1}))
        end
        args = args(l+1:end);
        if isempty(args), break, end
      end
      
      
      for itr = 1:obj.ntr
        if numel(color) == 3
          tmpColor = color;
          doColorbar = 0;
        elseif numel(color) == obj.ntr
          crange = [min(color) max(color)];
          if crange(1) == crange(2)
            crange = crange+crange(1)*0.1*[-1 1];
          end
          tmpColor = cmap2color(cmap,crange,color(itr));
          doColorbar = 1;
        elseif all(size(color) == [obj.ntr 3])
          tmpColor = color(itr,:);
          doColorbar = 1;
        end
        plot(h,obj(itr).x,obj(itr).z,'color',tmpColor);
        if doColorbar
          hcb = colorbar('peer',h);
          colormap(h,cmap)
          h.CLim = crange;
        end
        if itr == 1
          hold(gca,'on')
        end
      end
      hold(gca,'off')
      
      h.XGrid = 'on';
      h.YGrid = 'on';
      h.XLabel.String = 'x (d_i)';
      h.YLabel.String = 'z (d_i)';
    end
    function h = plot_all_xy(obj,varargin)
      % PICTRAJ.PLOT_ALL_XY Plots all trajectories in xy plane.
      h = obj.plot_all('x','y',varargin{:});
    end
    function h = plot_all_yz(obj,varargin)
      % PICTRAJ.PLOT_ALL_YZ Plots all trajectories in xy plane.      
      h = obj.plot_all('y','z',varargin{:});
    end
    function h = plot_all_zy(obj,varargin)
      % PICTRAJ.PLOT_ALL_YZ Plots all trajectories in xy plane.      
      h = obj.plot_all('z','y',varargin{:});
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