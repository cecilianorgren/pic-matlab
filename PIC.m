 classdef PIC
  % Load PIC simulation data
  %   Does not contain all the data, but loads it in an easily accesible manner  
  %
  %   pic = PIC(h5FilePath)
  %   Bx = pic.Bx; % Bx is a (nt x nx x ny) matrix
  %   B = pic.B; % structure with 3 (nt x nx x ny) matrices  
  
  properties (Access = protected)
    % Access = protected – access from class or subclasses
    % Data can be arbitrary size, so the class contains a pointer to the 
    % data file and each time loads the data with
    file_    
    info_
    fields_
    iteration_
    twpe_
    twci_
    xe_
    ze_
    xi_
    zi_
    grid_
    indices_
    it_
    ix_
    iz_
%    wpewce_ = [];
%    mime_ = [];
    
  end
  
  properties (Dependent = true)
    % Can be checked when setting, for example, right size, right type
    % Dependent properties don't store a value and can't be assigned
    % a value in their set method.
    file
    info
    fields
    iteration
    twpe
    twci
    xe
    ze
    xi
    zi
    grid
    indices
    it
    ix
    iz
    %wpewce
    %mime
    
  end
  
  properties (Constant = true)    
  end
  
  properties (Constant = true, Hidden = true)
    %MAX_TENSOR_ORDER = 2;
    %BASIS = {'xyz','xzy'}; % in order to use, need to define transformations between these
    %BASIS_NAMES = {'smilei','michael'};
  end

  properties (SetAccess = protected)    
    wpewce
    mime
    teti
    mass
  end
  
  properties
    userData = []; % anything can be added here
  end
  
  methods
    function obj = PIC(h5filePath)
      % sm = SMILEI(pathFields)
      % sm = SMILEI(pathFields,[],[]) - current implementation

      obj.file = h5filePath; 
      obj.info = h5info(h5filePath);
            
      obj.mass = obj.get_mass;
      uniqueMass = sort(unique(obj.mass));
      obj.mime = uniqueMass(2)/uniqueMass(1); % second lightest/lightest
      obj.teti = h5read(h5filePath,'/simulation_information/teti');
      obj.wpewce = h5read(h5filePath,'/simulation_information/wpewce');
      
      obj.iteration = get_iterations(obj);      
      obj.twpe = get_twpe(obj);      
      obj.twci = obj.twpe/(obj.wpewce*obj.mime);
      obj.indices_ = 1:numel(obj.iteration);
      
      obj.fields_ = get_fields(obj);
      obj.xe = h5read(h5filePath,'/simulation_information/xe'); % de
      obj.ze = h5read(h5filePath,'/simulation_information/ze'); % de
      obj.xi = obj.xe/sqrt(obj.mime);
      obj.zi = obj.ze/sqrt(obj.mime);
      obj.grid = {1:1:numel(obj.xe),1:1:numel(obj.ze)}; % originally, complete grid
      obj.ix = 1:1:numel(obj.xe);
      obj.iz = 1:1:numel(obj.ze);
      
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
          obj.iteration_ = builtin('subsref',obj.iteration,s);
          obj.twpe_ = builtin('subsref',obj.twpe,s);
          obj.twci_ = builtin('subsref',obj.twci,s);
          obj.indices_ = builtin('subsref',obj.indices,s);
          if numel(idx(1).subs) == 3 % time and two spatial indices
            s = substruct(idx(1).type,idx(1).subs(2));
            newgrid{1} = builtin('subsref',obj.grid{1},s); 
            obj.xe_ = builtin('subsref',obj.xe,s); 
            obj.xi_ = builtin('subsref',obj.xi,s); 
            s = substruct(idx(1).type,idx(1).subs(3));
            newgrid{2} = builtin('subsref',obj.grid{2},s);
            obj.ze_ = builtin('subsref',obj.ze,s); 
            obj.zi_ = builtin('subsref',obj.zi,s); 
            obj.grid_ = newgrid;
            %obj.gridsize_ = [numel(obj.grid{1}),numel(obj.grid{2})];      
          end
          
%           obj.iteration_ = builtin('subsref',obj.iteration,idx(1));
%           obj.twpe_ = builtin('subsref',obj.twpe,idx(1));
%           obj.twci_ = builtin('subsref',obj.twci,idx(1)); 
%           if numel(idx(1).subs)  1 % only time index
%             
%           end
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
    
    function value = length(obj)
      value = numel(obj.iteration);
    end
   
    % Get subset of data, time and/or space, using subsrefs for this for
    % now
    function obj = ilim(obj,value)
      % Get subset of output
      obj.twpe_ = obj.twpe_(value);
      obj.twci_ = obj.twci_(value);
      obj.iteration_ = obj.iteration_(value);      
      obj.indices_ = obj.indices_(value);      
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
    function obj = xlim(obj,value)
      % Get subset of x
      x1 = find(obj.xi >= value(1),1,'first'); 
      x2 = find(obj.xi <= value(2),1,'last'); 
      obj.xe_ = obj.xe_(x1:x2);
      obj.xi_ = obj.xi_(x1:x2);      
      obj.grid_{1} = obj.grid_{1}(x1:x2);
      obj.ix_ = obj.grid_{1};
    end
    function obj = zlim(obj,value)
      % Get subset of x
      z1 = find(obj.zi >= value(1),1,'first'); 
      z2 = find(obj.zi <= value(2),1,'last'); 
      obj.ze_ = obj.ze_(z1:z2);
      obj.zi_ = obj.zi_(z1:z2);      
      obj.grid_{2} = obj.grid_{2}(z1:z2);
      obj.iz_ = obj.grid_{2};
    end
    function obj = twpelim(obj,value)
      % Get subset of x
      i1 = find(obj.twpe >= value(1),1,'first'); 
      i2 = find(obj.twpe <= value(2),1,'last'); 
      obj.twpe_ = obj.twpe_(i1:i2);
      obj.twci_ = obj.twci_(i1:i2);
      obj.iteration_ = obj.iteration_(i1:i2);
    end
    function obj = twcilim(obj,value)
      % Get subset of x
      i1 = find(obj.twci >= value(1),1,'first'); 
      i2 = find(obj.twci <= value(2),1,'last'); 
      obj.twpe_ = obj.twpe_(i1:i2);
      obj.twci_ = obj.twci_(i1:i2);
      obj.iteration_ = obj.iteration_(i1:i2);
    end
    
    % Plotting routines, for simple diagnostics etc
    function [all_im, map] = make_gif(obj,fields,nrows,ncols,varargin)
      % [all_im, map] = MAKE_GIF(obj,fields,nrows,ncols)      
      % make gif
      % imwrite(im,map,'delme.gif','DelayTime',0.0,'LoopCount',0)  
      
      % Subsref error: Does not accept two outputs
      
      % Default options, values
      doAdjustCLim = 0;
      cmap = pic_colors('blue_red');
      doA = 0;
      
      nfields = numel(fields);
      ntimes = obj.length;
      
      nargs = numel(varargin);      
      if nargs > 0, have_options = 1; args = varargin(:); end
      
      while have_options
        l = 1;
        switch(lower(args{1}))
          case 'a'            
            if numel(args{2}) == 1              
              doA = args{2};
              levA = -25:1:0;
            else
              doA = 1;
              levA = args{2};
            end            
            args = args(l+1:end);
          case 'clim'
            l = 2;
            doAdjustCLim = 1;  
            clims = args{2};
            args = args(l+1:end);
          otherwise
            warning(sprintf('Input ''%s'' not recognized.',args{1}))
            args = args(l+1:end);
        end        
        if isempty(args), break, end    
      end
      
      % First load all data once and check what color limits we should use.
      %
      if 0% nfields == 1
        all_data = get_field(obj,fields{1});
        cmax = max(all_data(:));
        clim = cmax*[-1 1];
        doAdjustCLim = 1;
        cmap = pic_colors('blue_red');
      end
      
      % setup figure
      fig = figure;
      h = setup_subplots(nrows,ncols); % external function, must include in SMILEI.m
      disp('Adjust figure size, then hit any key to continue.')
      pause
      for itime = 1:ntimes
        for ifield = 1:nfields
          %tic;
          hca = h(ifield);
          S(1).type='()'; S(1).subs = {itime};
          data = get_field(obj.subsref(S),fields{ifield});
          imagesc(hca,obj.xi,obj.zi,squeeze(data)')
          hb = colorbar('peer',hca);
          hb.YLabel.String = fields{ifield};
          if doAdjustCLim
            hca.CLim = clims{ifield};            
            %colormap(cmap)
          end
          if doA
            hold(hca,'on')
            iAx = 1:4:obj.nx;
            iAz = 1:4:obj.nz;
            A = squeeze(get_field(obj.subsref(S),'A'));
            contour(hca,obj.xi(iAx),obj.zi(iAz),A(iAx,iAz)',levA,'k');
            hold(hca,'off')
          end
          hca.YDir = 'normal';
          hca.XLabel.String = 'x/d_i';
          hca.YLabel.String = 'z/d_i';
          %toc
        end
        if itime == 1
          colormap(cmap)
        end
        pause(0.1)
        if 1 % collect frames, for making gif
          iframe = itime;
          nframes = ntimes;
          currentBackgroundColor = get(gcf,'color');
          set(gcf,'color',[1 1 1]);
          drawnow      
          tmp_frame = getframe(gcf);
          %cell_movies{imovie}(itime) = tmp_frame;
          if iframe == 1 % initialize animated gif matrix
            [im_tmp,map] = rgb2ind(tmp_frame.cdata,256,'nodither');
            %map(end+1,:) = get(gcf,'color');
            im_tmp(1,1,1,nframes) = 0;                                                
            all_im = im_tmp;
          else
            all_im(:,:,1,iframe) = rgb2ind(tmp_frame.cdata,map,'nodither');
          end       
        end
      end
      %out = {all_im,map};
      
      % collect frames
      
    end
             
    
    % Data analysis routines, time derivatives etc.
    function out = calc_vector_potential(Bx,Bz,x,z)
      
    end
    
    % Get simulation meta data and parameters
    function out = get_twpe(obj)
      fileInfo = obj.info_;
      iGroup = find(contains({fileInfo.Groups.Name},'/data'));
      nIter = numel(fileInfo.Groups(iGroup).Groups); % number of iterations
      
      for iIter = 1:nIter
        % /data/00000xxxxx/ 
        % redo to actually find the 
        iAtt = find(contains({fileInfo.Groups(iGroup).Groups(iIter).Attributes.Name},'time'));
        time(iIter) = fileInfo.Groups(iGroup).Groups(iIter).Attributes(iAtt).Value;
      end
      out = time;
    end
    function out = get_iterations(obj)
      fileInfo = obj.info_;
      iGroup = find(contains({fileInfo.Groups.Name},'/data'));
      nOutput = numel(fileInfo.Groups(iGroup).Groups);
      for iOutput = 1:nOutput
        str = fileInfo.Groups(iGroup).Groups(iOutput).Name;
        split_str = strsplit(str,'/');
        iterations(iOutput) = str2num(split_str{3});
      end

      out = iterations;
    end
    function out = get_fields(obj)
      % needs to be adapted for the species subgroups
      fileInfo = obj.info;
      % fields structure is the same for all times
      out = {fileInfo.Groups(1).Groups(1).Datasets.Name};
      split_str = cellfun(@(x) strsplit(x,'/'),{fileInfo.Groups(1).Groups(1).Groups.Name},'UniformOutput',false);
      for iout = 1:numel(split_str)
        out{end+1} = split_str{iout}{end};
      end
      
    end
    function out = get_gridsize(obj)
      out = [numel(obj.xe) numel(obj.ze)];
    end
    function out = get_mass(obj)
      fileInfo = obj.info;
      iGroup = find(contains({fileInfo.Groups.Name},'/data'));
      iAtt = find(contains({fileInfo.Groups(iGroup).Groups(1).Groups(1).Datasets(1).Attributes.Name},'mass'));
      nSpecies = numel(fileInfo.Groups(iGroup).Groups(1).Groups(1).Datasets);
      for iSpecies = 1:nSpecies
        mass(iSpecies) = fileInfo.Groups(iGroup).Groups(1).Groups(1).Datasets(iSpecies).Attributes(iAtt).Value;
      end
      out = mass;
    end
    function out = get_charge(obj)
      fileInfo = obj.info;
      iGroup = find(contains({fileInfo.Groups.Name},'/data'));
      iAtt = find(contains({fileInfo.Groups(iGroup).Groups(1).Groups(1).Datasets(1).Attributes.Name},'charge'));
      nSpecies = numel(fileInfo.Groups(iGroup).Groups(1).Groups(1).Datasets);
      for iSpecies = 1:nSpecies
        charge(iSpecies) = fileInfo.Groups(iGroup).Groups(1).Groups(1).Datasets(iSpecies).Attributes(iAtt).Value;
      end
      out = charge;
    end
    function out = get_dfac(obj)
      fileInfo = obj.info;
      iGroup = find(contains({fileInfo.Groups.Name},'/data'));
      iAtt = find(contains({fileInfo.Groups(iGroup).Groups(1).Groups(1).Datasets(1).Attributes.Name},'dfac'));
      nSpecies = numel(fileInfo.Groups(iGroup).Groups(1).Groups(1).Datasets);
      for iSpecies = 1:nSpecies
        dfac(iSpecies) = fileInfo.Groups(iGroup).Groups(1).Groups(1).Datasets(iSpecies).Attributes(iAtt).Value;
      end
      out = dfac;
    end
    function out = nx(obj)
      out = numel(obj.xi);
    end
    function out = nz(obj)
      out = numel(obj.zi);
    end
    function out = nt(obj)
      out = obj.length;
    end
    
    % Get fields
    function out = A(obj)
      out = get_field(obj,'A');
    end
    function out = Bx(obj)
      out = get_field(obj,'bx')*obj.wpewce;
    end
    function out = By(obj)
      out = get_field(obj,'by')*obj.wpewce;
    end
    function out = Bz(obj)
      out = get_field(obj,'bz')*obj.wpewce;
    end
    function out = Ex(obj)
      out = get_field(obj,'ex')*sqrt(obj.mime)*obj.wpewce^2;
    end
    function out = Ey(obj)
      out = get_field(obj,'ey')*sqrt(obj.mime)*obj.wpewce^2;
    end
    function out = Ez(obj)
      out = get_field(obj,'ez')*sqrt(obj.mime)*obj.wpewce^2;
    end
    
    function out = ne(obj)
      % Get total electron density
      iSpecies = find(obj.get_charge == -1); % negatively charge particles are electrons
      dfac = obj.get_dfac;      
      var = zeros([1,obj.get_gridsize]);
      for iComp = 1:numel(iSpecies)
        dataset = sprintf('dns/%.0f',iSpecies(iComp));
        n_tmp = get_field(obj,dataset);
        var = var + n_tmp*dfac(iSpecies(iComp));
      end
      out = var;
    end
    function out = ni(obj)
      % Get total ion density
      iSpecies = find(obj.get_charge == 1); % negatively charge particles are electrons
      dfac = obj.get_dfac;      
      var = zeros([1,obj.get_gridsize]);
      for iComp = 1:numel(iSpecies)
        dataset = sprintf('dns/%.0f',iSpecies(iComp));
        n_tmp = get_field(obj,dataset);
        var = var + n_tmp*dfac(iSpecies(iComp));
      end
      out = var;
    end
    function out = n(obj,value)
      % Get total density of select species
      %   out = n(obj,value)
      dfac = obj.get_dfac;         
      dataset = sprintf('dns/%.0f',value);
      out = get_field(obj,dataset)*dfac(value);      
    end
    function out = jex(obj)
      % Get electron flux, x
      iSpecies = find(obj.get_charge == -1); % negatively charge particles are electrons
      dfac = obj.get_dfac;      
      var = zeros([1,obj.get_gridsize]);
      for iComp = 1:numel(iSpecies)
        dataset = sprintf('vxs/%.0f',iSpecies(iComp));
        n_tmp = get_field(obj,dataset);
        var = var + n_tmp*dfac(iComp)*obj.wpewce*sqrt(obj.mime);
      end
      out = var;
    end
    function out = jey(obj)
      % Get electron flux, y
      iSpecies = find(obj.get_charge == -1); % negatively charge particles are electrons
      dfac = obj.get_dfac;      
      var = zeros([1,obj.get_gridsize]);
      for iComp = 1:numel(iSpecies)
        dataset = sprintf('vys/%.0f',iSpecies(iComp));
        n_tmp = get_field(obj,dataset);
        var = var + n_tmp*dfac(iComp)*obj.wpewce*sqrt(obj.mime);
      end
      out = var;
    end
    function out = jez(obj)
      % Get electron flux, z
      iSpecies = find(obj.get_charge == -1); % negatively charge particles are electrons
      dfac = obj.get_dfac;      
      var = zeros([1,obj.get_gridsize]);
      for iComp = 1:numel(iSpecies)
        dataset = sprintf('vzs/%.0f',iSpecies(iComp));
        n_tmp = get_field(obj,dataset);
        var = var + n_tmp*dfac(iComp)*obj.wpewce*sqrt(obj.mime);
      end
      out = var;
    end
    function out = jix(obj)
      iSpecies = find(obj.get_charge == 1); % negatively charge particles are electrons
      dfac = obj.get_dfac;      
      var = zeros([1,obj.get_gridsize]);
      for iComp = 1:numel(iSpecies)
        dataset = sprintf('vxs/%.0f',iSpecies(iComp));
        n_tmp = get_field(obj,dataset);
        var = var + n_tmp*dfac(iComp)*obj.wpewce*sqrt(obj.mime);
      end
      out = var;
    end
    function out = jiy(obj)
      iSpecies = find(obj.get_charge == 1); % negatively charge particles are electrons
      dfac = obj.get_dfac;      
      var = zeros([1,obj.get_gridsize]);
      for iComp = 1:numel(iSpecies)
        dataset = sprintf('vys/%.0f',iSpecies(iComp));
        n_tmp = get_field(obj,dataset);
        var = var + n_tmp*dfac(iComp)*obj.wpewce*sqrt(obj.mime);
      end
      out = var;
    end
    function out = jiz(obj)
      iSpecies = find(obj.get_charge == 1); % negatively charge particles are electrons
      dfac = obj.get_dfac;      
      var = zeros([1,obj.get_gridsize]);
      for iComp = 1:numel(iSpecies)
        dataset = sprintf('vzs/%.0f',iSpecies(iComp));
        n_tmp = get_field(obj,dataset);
        var = var + n_tmp*dfac(iComp)*obj.wpewce*sqrt(obj.mime);
      end
      out = var;
    end
    function out = vex(obj)
      % Get electron velocity, x
      out = obj.jex./obj.ne;
    end
    function out = vey(obj)
      % Get electron velocity, y
      out = obj.jey./obj.ne;
    end
    function out = vez(obj)
      % Get electron velocity, z
      out = obj.jez./obj.ne;
    end
    function out = Jx(obj)
      out = obj.jix - obj.jex;
    end
    function out = Jy(obj)
      out = obj.jiy - obj.jey;
    end
    function out = Jz(obj)
      out = obj.jiz - obj.jez;
    end
    function out = vexx(obj)
      iSpecies = find(obj.get_charge == -1); % negatively charge particles are electrons
      dfac = obj.get_dfac;      
      var = zeros([1,obj.get_gridsize]);
      for iComp = 1:numel(iSpecies)
        dataset = sprintf('vxx/%.0f',iSpecies(iComp));
        var_tmp = get_field(obj,dataset);
        var = var + var_tmp*dfac(iComp)*mass(iSpecies)*wpewce^2;
      end
      out = var;
      out = [];
    end
    function out = vxx(obj,value)
      % Get total density of select species
      %   out = n(obj,value)
      dfac = obj.get_dfac;         
      dataset = sprintf('vxx/%.0f',value);
      out = obj.mass(value)*obj.wpewce^2*get_field(obj,dataset)*dfac(value);      
    end
    function out = vyy(obj,value)
      % Get total density of select species
      %   out = n(obj,value)
      dfac = obj.get_dfac;         
      dataset = sprintf('vyy/%.0f',value);
      out = obj.mass(value)*obj.wpewce^2*get_field(obj,dataset)*dfac(value);      
    end
    function out = vzz(obj,value)
      % Get total density of select species
      %   out = n(obj,value)
      dfac = obj.get_dfac;         
      dataset = sprintf('vzz/%.0f',value);
      out = obj.mass(value)*obj.wpewce^2*get_field(obj,dataset)*dfac(value);      
    end
    function out = vv_diag(obj,value)
      % Get 
      %   out = vv_diag(obj,value)
      dfac = obj.get_dfac;         
      vxx = get_field(obj,sprintf('vxx/%.0f',value))*dfac(value);
      vyy = get_field(obj,sprintf('vyy/%.0f',value))*dfac(value);
      vzz = get_field(obj,sprintf('vzz/%.0f',value))*dfac(value);
      out = obj.mass(value)*obj.wpewce^2*(vxx + vyy + vzz)/3;
    end
    function out = PB(obj)
      % Get total density of select species
      %   out = n(obj,value)      
      Bx = obj.Bx;
      By = obj.By;
      Bz = obj.Bz;
      out = 0.5*sqrt(Bx.^2 + By.^2 + Bz.^2);      
    end
    
    function [vxx,vxy,vxz,vyy,vyz,vzz] = vv(obj,iSpecies_orig)
      % [vxx,vxy,vxz,vyy,vyz,vzz] = vv(obj,iSpecies_orig)
      %nargout      
      mass_sp = obj.mass; 
      if numel(unique(mass_sp(iSpecies_orig))) > 1
        error('All species do not have the same mass.'); 
      else
        mass_sp = mass_sp(iSpecies_orig(1));
      end        
      dfac = obj.get_dfac;
      nSpecies = numel(iSpecies_orig);
      
      vxx = zeros(obj.length,obj.nx,obj.nz);
      vyy = zeros(obj.length,obj.nx,obj.nz);
      vzz = zeros(obj.length,obj.nx,obj.nz);
      vxy = zeros(obj.length,obj.nx,obj.nz);
      vxz = zeros(obj.length,obj.nx,obj.nz);
      vyz = zeros(obj.length,obj.nx,obj.nz);
        
      for iSpecies = iSpecies_orig
        dset_vxx = sprintf('vxx/%.0f',iSpecies);
        dset_vxy = sprintf('vxy/%.0f',iSpecies);
        dset_vxz = sprintf('vxz/%.0f',iSpecies);
        dset_vyy = sprintf('vyy/%.0f',iSpecies);
        dset_vyz = sprintf('vyz/%.0f',iSpecies);
        dset_vzz = sprintf('vzz/%.0f',iSpecies);
        vxx = vxx + obj.get_field(dset_vxx)*dfac(iSpecies)*mass_sp*obj.wpewce^2;
        vxy = vxy + obj.get_field(dset_vxy)*dfac(iSpecies)*mass_sp*obj.wpewce^2;
        vxz = vxz + obj.get_field(dset_vxz)*dfac(iSpecies)*mass_sp*obj.wpewce^2;
        vyy = vyy + obj.get_field(dset_vyy)*dfac(iSpecies)*mass_sp*obj.wpewce^2;
        vyz = vyz + obj.get_field(dset_vyz)*dfac(iSpecies)*mass_sp*obj.wpewce^2;
        vzz = vzz + obj.get_field(dset_vzz)*dfac(iSpecies)*mass_sp*obj.wpewce^2;         
      end
      
      
    end
    function [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = njp(obj,iSpecies)
      % [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = njp(obj,iSpecies)
      iSpecies_orig = iSpecies;
      nSpecies = numel(iSpecies);
      dfac = obj.get_dfac;
      
      % Get density, flux and pressure of given species
      if nSpecies == 1 % only one species
        % n
        dset_n = sprintf('dns/%.0f',iSpecies);
        n = obj.get_field(dset_n)*dfac(iSpecies);
        % j
        dset_jx = sprintf('vxs/%.0f',iSpecies);
        dset_jy = sprintf('vys/%.0f',iSpecies);
        dset_jz = sprintf('vzs/%.0f',iSpecies);
        vxs = obj.get_field(dset_jx)*dfac(iSpecies);
        vys = obj.get_field(dset_jy)*dfac(iSpecies);
        vzs = obj.get_field(dset_jz)*dfac(iSpecies);      
        jx = vxs*obj.wpewce*sqrt(obj.mime);
        jy = vys*obj.wpewce*sqrt(obj.mime);
        jz = vzs*obj.wpewce*sqrt(obj.mime);
        % p
        dset_vxx = sprintf('vxx/%.0f',iSpecies);
        dset_vxy = sprintf('vxy/%.0f',iSpecies);
        dset_vxz = sprintf('vxz/%.0f',iSpecies);
        dset_vyy = sprintf('vyy/%.0f',iSpecies);
        dset_vyz = sprintf('vyz/%.0f',iSpecies);
        dset_vzz = sprintf('vzz/%.0f',iSpecies);
        vxx = obj.get_field(dset_vxx)*dfac(iSpecies);
        vxy = obj.get_field(dset_vxy)*dfac(iSpecies);
        vxz = obj.get_field(dset_vxz)*dfac(iSpecies);
        vyy = obj.get_field(dset_vyy)*dfac(iSpecies);
        vyz = obj.get_field(dset_vyz)*dfac(iSpecies);
        vzz = obj.get_field(dset_vzz)*dfac(iSpecies);       
          
        pxx = obj.mass(iSpecies)*obj.wpewce^2*( vxx - vxs.*vxs./n );
        pxy = obj.mass(iSpecies)*obj.wpewce^2*( vxy - vxs.*vys./n );
        pxz = obj.mass(iSpecies)*obj.wpewce^2*( vxz - vxs.*vzs./n );
        pyy = obj.mass(iSpecies)*obj.wpewce^2*( vyy - vys.*vys./n );
        pyz = obj.mass(iSpecies)*obj.wpewce^2*( vyz - vys.*vzs./n );
        pzz = obj.mass(iSpecies)*obj.wpewce^2*( vzz - vzs.*vzs./n );
      else % group multiple species
        % check so that all species have the same mass and charge
        mass_sp = obj.mass; 
        if numel(unique(mass_sp(iSpecies_orig))) > 1
          error('All species do not have the same mass.'); 
        else
          mass_sp = mass_sp(iSpecies_orig(1));
        end
        charge = obj.get_charge; charge(iSpecies_orig);
        if numel(unique(charge(iSpecies_orig))) > 1
          error('All species do not have the same charge.');         
        end
        
        n = zeros(obj.length,numel(obj.xi),numel(obj.zi));
        jx = zeros(obj.length,numel(obj.xi),numel(obj.zi));
        jy = zeros(obj.length,numel(obj.xi),numel(obj.zi));
        jz = zeros(obj.length,numel(obj.xi),numel(obj.zi));
        vxs = zeros(obj.length,numel(obj.xi),numel(obj.zi));
        vys = zeros(obj.length,numel(obj.xi),numel(obj.zi));
        vzs = zeros(obj.length,numel(obj.xi),numel(obj.zi));
        pxx = zeros(obj.length,numel(obj.xi),numel(obj.zi));
        pyy = zeros(obj.length,numel(obj.xi),numel(obj.zi));
        pzz = zeros(obj.length,numel(obj.xi),numel(obj.zi));
        pxy = zeros(obj.length,numel(obj.xi),numel(obj.zi));
        pxz = zeros(obj.length,numel(obj.xi),numel(obj.zi));
        pyz = zeros(obj.length,numel(obj.xi),numel(obj.zi));
        vxx = zeros(obj.length,numel(obj.xi),numel(obj.zi));
        vyy = zeros(obj.length,numel(obj.xi),numel(obj.zi));
        vzz = zeros(obj.length,numel(obj.xi),numel(obj.zi));
        vxy = zeros(obj.length,numel(obj.xi),numel(obj.zi));
        vxz = zeros(obj.length,numel(obj.xi),numel(obj.zi));
        vyz = zeros(obj.length,numel(obj.xi),numel(obj.zi));
        
        for iSpecies = iSpecies_orig
          % n 
          % number of macroparticles, this can be added straight up
          dset_n = sprintf('dns/%.0f',iSpecies);
          n = n + obj.get_field(dset_n)*dfac(iSpecies);
          % v (j)
          dset_jx = sprintf('vxs/%.0f',iSpecies);
          dset_jy = sprintf('vys/%.0f',iSpecies);
          dset_jz = sprintf('vzs/%.0f',iSpecies);
          vxs = vxs + obj.get_field(dset_jx)*dfac(iSpecies);
          vys = vys + obj.get_field(dset_jy)*dfac(iSpecies);
          vzs = vzs + obj.get_field(dset_jz)*dfac(iSpecies);          
          % vv
          dset_vxx = sprintf('vxx/%.0f',iSpecies);
          dset_vxy = sprintf('vxy/%.0f',iSpecies);
          dset_vxz = sprintf('vxz/%.0f',iSpecies);
          dset_vyy = sprintf('vyy/%.0f',iSpecies);
          dset_vyz = sprintf('vyz/%.0f',iSpecies);
          dset_vzz = sprintf('vzz/%.0f',iSpecies);
          vxx = vxx + obj.get_field(dset_vxx)*dfac(iSpecies);
          vxy = vxy + obj.get_field(dset_vxy)*dfac(iSpecies);
          vxz = vxz + obj.get_field(dset_vxz)*dfac(iSpecies);
          vyy = vyy + obj.get_field(dset_vyy)*dfac(iSpecies);
          vyz = vyz + obj.get_field(dset_vyz)*dfac(iSpecies);
          vzz = vzz + obj.get_field(dset_vzz)*dfac(iSpecies);
        end
        % j = nv
        jx = vxs*obj.wpewce*sqrt(obj.mime);
        jy = vys*obj.wpewce*sqrt(obj.mime);
        jz = vzs*obj.wpewce*sqrt(obj.mime);
        % p = P - mnvv
        pxx = mass_sp*obj.wpewce^2*(vxx - vxs.*vxs./n); 
        pxy = mass_sp*obj.wpewce^2*(vxy - vxs.*vys./n);
        pxz = mass_sp*obj.wpewce^2*(vxz - vxs.*vzs./n);
        pyy = mass_sp*obj.wpewce^2*(vyy - vys.*vys./n);
        pyz = mass_sp*obj.wpewce^2*(vyz - vys.*vzs./n);
        pzz = mass_sp*obj.wpewce^2*(vzz - vzs.*vzs./n);
      end
    end
    
    % Ge derived quantities
    function out = UB(obj)
      % Magnetic energy density 0.5*(Bx^2 + By^2 + Bz^2) summed up 
      out = h5read(obj.file,'/scalar_timeseries/U/B');
      out = out(obj.indices_);
    end
    function out = dUB(obj)
      % Magnetic energy density 0.5*(Bx^2 + By^2 + Bz^2) summed up 
      out = h5read(obj.file,'/scalar_timeseries/U/B');
      out = abs([0 cumsum(diff(out-out(1)))]);
      out = out(obj.indices_);
    end
    function out = UK(obj,value)
      % 
      out = h5read(obj.file,['/scalar_timeseries/U/K/' num2str(value)]);
      out = out(obj.indices_);
    end
    function out = UT(obj,value)
      % 
      out = h5read(obj.file,['/scalar_timeseries/U/T/' num2str(value)]);
      out = out(obj.indices_);
    end
    function out = RE(obj)
      % Reconnection rate from out-of-plane electric field Ey at X line
      out = h5read(obj.file,'/scalar_timeseries/R/Ey');
      out = out(obj.indices_);
    end
    function out = RA(obj)
      % Reconnection rate from vector potential dA/dt at X line
      out = h5read(obj.file,'/scalar_timeseries/R/A');
      out = out(obj.indices_);
    end        
    
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
    function obj = set.iteration(obj,value)
      obj.iteration_ = value;
    end
    function obj = set.twpe(obj,value)
      obj.twpe_ = value;
    end
    function obj = set.twci(obj,value)
      obj.twci_ = value;
    end
    function obj = set.xe(obj,value)
      obj.xe_ = value;
    end
    function obj = set.ze(obj,value)
      obj.ze_ = value;
    end
    function obj = set.xi(obj,value)
      obj.xi_ = value;
    end
    function obj = set.zi(obj,value)
      obj.zi_ = value;
    end
    function obj = set.grid(obj,value)
      obj.grid_ = value;
    end
    function obj = set.indices(obj,value)
      obj.indices_ = value;
    end
    function obj = set.it(obj,value)
      obj.it_ = value;
    end
    function obj = set.ix(obj,value)
      obj.ix_ = value;
    end
    function obj = set.iz(obj,value)
      obj.iz_ = value;
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
    function value = get.iteration(obj)
      value = obj.iteration_;
    end 
    function value = get.twpe(obj)
      value = obj.twpe_;
    end
    function value = get.twci(obj)
      value = obj.twci_;
    end 
    function value = get.xe(obj)
      value = obj.xe_;
    end
    function value = get.ze(obj)
      value = obj.ze_;
    end
    function value = get.xi(obj)
      value = obj.xi_;
    end
    function value = get.zi(obj)
      value = obj.zi_;
    end
    function value = get.grid(obj)
      value = obj.grid_;
    end
    function value = get.indices(obj)
      value = obj.indices_;
    end
    function value = get.it(obj)
      value = obj.it_;
    end
    function value = get.ix(obj)
      value = obj.ix_;
    end
    function value = get.iz(obj)
      value = obj.iz_;
    end
  end
  
  methods (Access = protected)
    function out = get_field(obj,field)
      % get iterations
      iterations = obj.iteration;
      nIter = obj.length;
      % initialize matrix
      data = nan([nIter,obj.get_gridsize]);
      for iIter = 1:nIter
        iter = iterations(iIter);
        str_iter = sprintf('%010.0f',iter);
        if 1 % read part of data, seems to be slightly faster
          data_tmp = h5read(obj.file,...
             ['/data/' str_iter '/' field],...
             [obj.grid{1}(1) obj.grid{2}(1)],... % start indices
             [numel(obj.grid{1}) numel(obj.grid{2})]); % number of counts
          data(iIter,:,:) = data_tmp;%(obj.grid{1},obj.grid{2});
        else % read all data then pick indices
          data_tmp = h5read(obj.file,['/data/' str_iter '/' field]); % de;
          data(iIter,:,:) = data_tmp(obj.grid{1},obj.grid{2});        
        end
      end
      out = data;
    end
    
    function Ts = changeBasis(obj, flag)
      % Tranform from one coordinate system to another and return new
      % TimeSeries.
      % flag: = 'xyz>rlp' - Cartesian XYZ to spherical latitude
      %         'rlp>xyz' - Spherical latitude to cartesian XYZ
      %         'xyz>rpz' - Cartesian XYZ to cylindrical
      %         'rpz>xyz' - Cylidrical to cartesian XYZ
      %         'xyz>rtp' - Cartesian XYZ to spherical colatitude
      %         'rtp>xyz' - Spherical colatitude to cartesian XYZ
      %         'rtp>rlp' - Spherical colatitude to spherical latitude
      %         'rlp>rtp' - Spherical latitude to colatitude
      switch lower(flag)
        case 'xyz>rlp'
          [phi, lambda, r] = cart2sph(obj.x.data, obj.y.data, obj.z.data);
          Ts = TSeries(obj.time, [r, lambda, phi], 'vec_rlp');
        case 'rlp>xyz'
          [x, y, z] = sph2cart(obj.phi.data, obj.lambda.data, obj.r.data);
          Ts = TSeries(obj.time, [x, y, z], 'vec_xyz');
        case 'xyz>rpz'
          [phi, r, z] = cart2pol(obj.x.data, obj.y.data, obj.z.data);
          Ts = TSeries(obj.time, [r, phi, z], 'vec_rpz');
        case 'rpz>xyz'
          [x, y, z] = pol2cart(obj.phi.data, obj.r.data, obj.z.data);
          Ts = TSeries(obj.time, [x, y, z], 'vec_xyz');
        case 'xyz>rtp'
          [phi, lambda, r] = cart2sph(obj.x.data, obj.y.data, obj.z.data);
          theta = pi/2 - lambda;
          Ts = TSeries(obj.time, [r, theta, phi], 'vec_rtp');
        case 'rtp>xyz'
          lambda = pi/2 - obj.theta.data;
          [x, y, z] = sph2cart(obj.phi.data, lambda, obj.r.data);
          Ts = TSeries(obj.time, [x, y, z], 'vec_xyz');
        case 'rtp>rlp'
          lambda = pi/2 - obj.theta.data;
          Ts = TSeries(obj.time, [obj.r.data,lambda,obj.phi.data],'vec_rlp');
        case 'rlp>rtp'
          theta = pi/2 - obj.lambda.data;
          Ts = TSeries(obj.time, [obj.r.data,theta,obj.phi.data],'vec_rtp');
        case 'xy>rp'
          [phi, r] = cart2pol(obj.x.data, obj.y.data);
          Ts = TSeries(obj.time, [r, phi], 'vec_rp');
        case 'rp>xy'
          [x, y] = pol2cart(obj.phi.data, obj.r.data);
          Ts = TSeries(obj.time, [x, y], 'vec_xy');
        otherwise
          errStr='Invalid transformation'; error(errStr);
      end
    end
  end
  
end