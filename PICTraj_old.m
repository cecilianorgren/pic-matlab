classdef PICTraj_
  % Load PIC trajectories data
  %   Does not contain all the data, but loads it in an easily accesible manner  
  %
  %   tr = PICTraj(h5FilePath)
  
  properties (Access = protected)
    % Access = protected – access from class or subclasses
    % Data can be arbitrary size, so the class contains a pointer to the 
    % data file and each time loads the data with
    file_
    info_
    fields_
    attributes_
    id_
    twpe_
    twci_
    xi_
    zi_
    it_
    nt_
  end
  
  properties (Dependent = true)
    % Can be checked when setting, for example, right size, right type
    % Dependent properties don't store a value and can't be assigned
    % a value in their set method.
    file
    info
    fields
    attributes
    trajectories % indices od trajectories
    id
    twpe
    twci
    xi
    zi  
    it
    nt
  end
  
  methods
    function obj = PICTraj(h5filePath)
      % traj = PICTraj(pathFields)
            
      obj.file = h5filePath; 
      obj.info = h5info(h5filePath);   
      obj.nt_ = obj.get_ntraj;
      obj.fields_ = obj.get_fields;
      obj.attributes_ = obj.get_attributes;
      obj.it = 1:obj.nt_;
    end
    
    function [varargout] = subsref(obj,idx)
      %SUBSREF handle indexing
%      nargout
%      idx(:)
      switch idx(1).type
        % Use the built-in subsref for dot notation
        case '.'
          [varargout{1:nargout}] = builtin('subsref',obj,idx);
        case '()'
          %nargout
          % first index is time
          s = substruct(idx(1).type,idx(1).subs(1));
          obj.it_ = builtin('subsref',obj.it_,s);
          obj.nt_ = numel(obj.it_);
          
          if numel(idx) > 1
            obj = builtin('subsref',obj,idx(2:end));
          end
          try
          [varargout{1:nargout}] = obj;
          catch
            1;
          end
        case '{}'
          error('SMILEI:subsref',...
            'Not a supported subscripted reference.')
      end
    end  
    
    function value = ntr(obj)
      value = numel(obj.id);
    end
   
    % Get subset of data, time and/or space, using subsrefs for this for
    % now
    function obj = ilim(obj,value)
      % Get subset of output
      %obj.twpe_ = obj.twpe_(value);
      %obj.twci_ = obj.twci_(value);
      obj.iteration_ = obj.iteration_(value);      
      %obj.indices_ = obj.indices_(value);      
    end
    function obj = i(obj,varargin)
      % Get subset of output
      nargs = numel(varargin);
      switch nargs 
        case 1 % time
          it = varargin{1};
          obj = obj.ilim(it);
        case 3 % time, x, z
          it = varargin{1};
          ix = varargin{2};
          iz = varargin{3};
          obj = obj.ilim(it);
          obj.xe_ = obj.xe_(ix);
          obj.xi_ = obj.xi_(ix);      
          obj.grid_{1} = obj.grid_{1}(ix);
          obj.ze_ = obj.ze_(iz);
          obj.zi_ = obj.zi_(iz);      
          obj.grid_{2} = obj.grid_{2}(iz);
        otherwise
          error('Number of inputs not supported.')
      end          
    end
    function obj = xlim(obj,value,varargin)
      % Get subset of x
      nt = obj.nt;
      nd = obj.nd;
      method = 'center';
      
      for it = 1:nt
        xtmp = obj.xi{it};
        xtmp = [xtmp{:}];
        xtmp = reshape(xtmp,[2,nd{it}])'; % nx x nd matrix
        if strcmp(method,'center')
          xtmp = (xtmp(:,1)+xtmp(:,2))/2*[1 1];
        end
        i1 = find(xtmp(:,1) >= value(1)); 
        i2 = find(xtmp(:,2) <= value(2)); 
        inds{it} = intersect(i1,i2);
      end
            
      %nx = cellfun(@(x) numel(x),obj.dists,'UniformOutput',false);      
      %inds = obj.ind_from_lim(obj.xi_,value,varargin{:});
      
      % Update grid and indices
      obj = obj.update_inds(inds);      
    end
    function obj = zlim(obj,value,varargin)
      % Get subset of z
      nt = obj.nt;
      nd = obj.nd;
      method = 'center';
      %method = 'edges';
      
      for it = 1:nt
        xtmp = obj.zi{it};
        xtmp = [xtmp{:}];
        xtmp = reshape(xtmp,[2,nd{it}])'; % nx x nd matrix
        if strcmp(method,'center')
          xtmp = (xtmp(:,1)+xtmp(:,2))/2*[1 1];
        end
        i1 = find(xtmp(:,1) >= value(1)); 
        i2 = find(xtmp(:,2) <= value(2)); 
        inds{it} = intersect(i1,i2);
      end
      
      % Update grid and indices
      obj = obj.update_inds(inds);   
    end
    function obj = twpelim(obj,value,varargin)
      % Get subset of twpe
      inds = obj.ind_from_lim(obj.twpe_,value,varargin{:});
      obj.twpe_ = obj.twpe_(inds);
      obj.twci_ = obj.twci_(inds);
      obj.it_ = obj.it_(inds);
      obj.iteration_ = obj.iteration_(inds);
    end
    function obj = twcilim(obj,value,varargin)
      % Get subset of twci
      inds = obj.ind_from_lim(obj.twci,value,varargin{:});
      obj.twpe_ = obj.twpe_(inds);
      obj.twci_ = obj.twci_(inds);
      obj.it_ = obj.it_(inds);
      obj.indices_ = obj.indices_(inds);
      obj.iteration_ = obj.iteration_(inds);
    end
    function obj = xpass(obj,value)
      % choose subset of data that pass.
    end
    
    function out = get_init(obj)
      % obj.get_init
      % get [x0,y0,z0,vx0,vy0,vz0,t0]
      
      init = zeros(obj.nt,7);
      for it = 1:obj.nt
        iTraj = obj.it(it);
        str_traj = sprintf('%06.0f',iTraj);
        init(it,1) = h5readatt(obj.file,['/traj/' str_traj],'x0');
        init(it,2) = h5readatt(obj.file,['/traj/' str_traj],'y0');
        init(it,3) = h5readatt(obj.file,['/traj/' str_traj],'z0');
        init(it,4) = h5readatt(obj.file,['/traj/' str_traj],'vx0');
        init(it,5) = h5readatt(obj.file,['/traj/' str_traj],'vy0');
        init(it,6) = h5readatt(obj.file,['/traj/' str_traj],'vz0');
        init(it,7) = h5readatt(obj.file,['/traj/' str_traj],'t0');
      end     
      out = init;
    end
    function out = t0(obj)
      % obj.x0
      % get x0
      out = obj.read_att('t0');      
    end
    function out = x0(obj)
      % obj.x0
      % get x0
      out = obj.read_att('x0');      
    end
    function out = y0(obj)
      % obj.x0
      % get x0
      out = obj.read_att('x0');      
    end
    function out = z0(obj)
      % obj.x0
      % get x0
      out = obj.read_att('y0');      
    end
    function out = vx0(obj)
      % obj.x0
      % get x0
      out = obj.read_att('vx0');      
    end
    function out = vy0(obj)
      % obj.x0
      % get x0
      out = obj.read_att('vy0');      
    end
    function out = vz0(obj)
      % obj.x0
      % get x0
      out = obj.read_att('vz0');      
    end
    function out = get_fields(obj)
      fileInfo = obj.info;
      % fields structure is the same for all times
      out = {fileInfo.Groups(1).Groups(1).Datasets.Name};      
    end
    function out = get_attributes(obj)
      fileInfo = obj.info;
      % fields structure is the same for all times
      out = {fileInfo.Groups(1).Groups(1).Attributes.Name};
    end
    function out = get_ntraj(obj)
      fileInfo = obj.info;
      out = numel(fileInfo.Groups(1).Groups);
    end
    % Plotting routines, for simple diagnostics etc
    
    % Load corresponding phase space distribution
    function varargout = fxyz(obj,it,id,iss,sumdim)
            
      iSpecies = iss;        
      dataset = obj.dataset_str(it,id);
      info = h5info(obj.file,dataset);      
      nSpecies = info.Dataspace.Size(4);
      datasize = info.Dataspace.Size(1:3);
            
      allAxes = h5readatt(obj.file,dataset,'axes');
      selAxes = allAxes(:,iSpecies);
      nUniqueAxes = size(unique(selAxes','rows'),1);
      comAxes = min(selAxes(:)):min(min(diff(selAxes,1))):max(selAxes(:)); % keep max resolution and max axes value
      newDataSize = [numel(comAxes) numel(comAxes) numel(comAxes)];
      
      dist = zeros(newDataSize);
      for iSpeciesTmp = iSpecies
        data_tmp = h5read(obj.file,dataset,[1 1 1,iSpeciesTmp],[datasize,1]);
        if nUniqueAxes > 1
          [Vx,Vy,Vz] = meshgrid(allAxes(:,iSpeciesTmp),allAxes(:,iSpeciesTmp),allAxes(:,iSpeciesTmp));
          [Vxq,Vyq,Vzq] = meshgrid(comAxes,comAxes,comAxes);
          data_tmp = interp3(Vx,Vy,Vz,data_tmp,Vxq,Vyq,Vzq);
          data_tmp(isnan(data_tmp)) = 0;
        end                    
        dist = dist + data_tmp;
      end
            
      % Reduce distribution
      if exist('sumdim','var') && any(sumdim == [1 2 3])
        for isum = numel(sumdim):-1:1
          dist = squeeze(sum(dist,sumdim(isum)));
        end
      end
      
      x = h5readatt(obj.file,dataset,'x');
      z = h5readatt(obj.file,dataset,'z');
      
      if nargout == 1
        out.f = dist;
        out.v = comAxes;
        out.x = x;
        out.z = z;
        varargout{1} = out;
      elseif nargout == 3
        varargout{1} = comAxes;
        varargout{2} = comAxes;
        varargout{1} = dist;
      end
    end
    function varargout = f(obj,it,id,iss)
            
      iSpecies = iss;        
      dataset = obj.dataset_str(it,id);
      info = h5info(obj.file,dataset);      
      nSpecies = info.Dataspace.Size(4);
      datasize = info.Dataspace.Size(1:3);
            
      allAxes = h5readatt(obj.file,dataset,'axes');
      selAxes = allAxes(:,iSpecies);
      nUniqueAxes = size(unique(selAxes','rows'),1);
      comAxes = min(selAxes(:)):min(min(diff(selAxes,1))):max(selAxes(:)); % keep max resolution and max axes value
      newDataSize = [numel(comAxes) numel(comAxes) numel(comAxes)];
      
      dist = zeros(newDataSize);
      for iSpeciesTmp = iSpecies
        data_tmp = h5read(obj.file,dataset,[1 1 1,iSpeciesTmp],[datasize,1]);
        if nUniqueAxes > 1
          [Vx,Vy,Vz] = meshgrid(allAxes(:,iSpeciesTmp),allAxes(:,iSpeciesTmp),allAxes(:,iSpeciesTmp));
          [Vxq,Vyq,Vzq] = meshgrid(comAxes,comAxes,comAxes);
          data_tmp = interp3(Vx,Vy,Vz,data_tmp,Vxq,Vyq,Vzq);
          data_tmp(isnan(data_tmp)) = 0;
        end                    
        dist = dist + data_tmp;
      end
            
      % Reduce distribution
      fstr = {'yz','xz','xy'};
      for sumdim = 1:3        
        for isum = numel(sumdim):-1:1
          dist_sum{sumdim} = squeeze(sum(dist,sumdim(isum)));
        end        
      end
      
      x = h5readatt(obj.file,dataset,'x');
      z = h5readatt(obj.file,dataset,'z');
      
      out.f = dist;
      out.fyz = dist_sum{1};
      out.fxz = dist_sum{2};
      out.fxy = dist_sum{3};
      out.v = comAxes;
      out.x = x;
      out.z = z;
      varargout{1} = out;
    end        
    
    % Plotting functions
    
    % Interpolate distribution
    function out = interp(obj,varstr,x,z,t,varargin)
      % Interpolates fields to given x,z,t
      
      % Check input
      
      % Check if inteprolation is needed
      
      % Load bounding data
   %   var = 
      
      % Interpolate
      switch method
        case 1
      end
    end
    
    % Ge derived quantities      
    
    % Get and set properties    
    function obj = set.info(obj,value)
      obj.info_ = value;
    end
    function obj = set.file(obj,value)
      obj.file_ = value;
    end
    function obj = set.fields(obj,value)
      obj.fields_ = value;
    end
    function obj = set.attributes(obj,value)
      obj.attributes_ = value;
    end
    function obj = set.twpe(obj,value)
      obj.twpe_ = value;
    end
    function obj = set.twci(obj,value)
      obj.twci_ = value;
    end
    function obj = set.xi(obj,value)
      obj.xi_ = value;
    end
    function obj = set.zi(obj,value)
      obj.zi_ = value;
    end
    function obj = set.it(obj,value)
      obj.it_ = value;
    end    
    function obj = set.nt(obj,value)
      obj.nt_ = value;
    end
    
    function value = get.info(obj)
      value = obj.info_;
    end
    function value = get.file(obj)
      value = obj.file_;
    end
    function value = get.fields(obj)
      value = obj.fields_;
    end  
    function value = get.attributes(obj)
      value = obj.attributes_;
    end  
    function value = get.twpe(obj)
      value = obj.twpe_;
    end
    function value = get.twci(obj)
      value = obj.twci_;
    end 
    function value = get.xi(obj)
      value = obj.xi_;
    end
    function value = get.zi(obj)
      value = obj.zi_;
    end
    function value = get.it(obj)
      value = obj.it_;
    end
    function value = get.nt(obj)
      value = obj.nt_;
    end   
  end
  
  methods (Static) % does not require object as input, but still needs to be called as obj.func (?)
    function out = ind_from_lim(var,value,varargin)
      % method is the same for xlim, zlim ilim, i, twpelim, twcilim
      
      % Defaults
      doBounding = 0;
      nBounding = 0;
      doClosest = 0;
      nClosest = 1; % only the closest index
      
      if numel(value) == 1
        doClosest = 1;        
      end
      
      % Check input
      have_options = 0;
      nargs = numel(varargin);      
      if nargs > 0, have_options = 1; args = varargin(:); end
      
      while have_options
        l = 1;
        switch(lower(args{1}))
          case 'closest'
            l = 2;
            doClosest = 1;
            nClosest = args{2};            
            args = args(l+1:end);    
          case 'bounding'
            l = 2;
            doBounding = 1;
            nBounding = args{2};            
            args = args(l+1:end);        
          otherwise
            warning(sprintf('Input ''%s'' not recognized.',args{1}))
            args = args(l+1:end);
        end        
        if isempty(args), break, end    
      end
      
      
      % Find indices
      if doBounding
        i0 = find(abs(var-value(1)) == min(abs(var-value(1))));
        i1 = i0 - nBounding;
        i2 = i0 + nBounding;
        % Check so that indices are not outside range
        if i1 < 1
          i1 = 1; 
        end
        if i2 > numel(var) 
          i2 = numel(var) ; 
        end
        inds = i1:i2;
      elseif doClosest
        ii = abs(var-value(1));
        [is, index] = sort(abs(var-value(1)));
        inds = sort(index(1:nClosest));
      else
        i1 = find(var >= value(1),1,'first'); 
        i2 = find(var <= value(2),1,'last'); 
        inds = i1:i2;
      end
      
      out = inds;
    end        
  end
  
  methods (Access = protected)
    function out = read_trajectories(obj)
      % get indices      
      % initialize empty structure
      data = struct;
      data_fields = cell(obj.nt,numel(obj.fields));
      for it = 1:obj.nt
        iTraj = obj.it(it);
        str_traj = sprintf('%06.0f',iTraj);
        for iField = 1:numel(obj.fields)
          data_fields{it,iField} = h5read(obj.file,['/traj/' str_traj '/' obj.fields{iField}]);
        end        
      end
      out = cell2struct(data_fields,obj.fields,2);
    end    
    function out = read_fields(obj,field)
      % get indices      
      % initialize empty structure
      data = struct;
      for iTraj = obj.it        
        str_traj = sprintf('%06.0f',iTraj);        
        data_tmp = h5read(obj.file,...
           ['/data/' str_traj '/' field]);        
        %eval('data.'
      end
      out = data;
    end
    
    function out = read_att(obj,att)
      % obj.get_init
      % get [x0,y0,z0,vx0,vy0,vz0,t0]
      
      init = zeros(obj.nt,1);
      for it = 1:obj.nt
        iTraj = obj.it(it);
        str_traj = sprintf('%06.0f',iTraj);
        init(it,1) = h5readatt(obj.file,['/traj/' str_traj],att);
      end     
      out = init;
    end
  end
  
end