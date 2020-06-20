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
      TR(itr)
    end
    % Finding subsets of data
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
    function h = plot_all_xz(obj,varargin)
      %hca = axes;
      
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
      %hca = axes;
      
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
      %hca = axes;
      
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