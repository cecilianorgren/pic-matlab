classdef PIC3D
  % Load PIC simulation data
  %   Does not contain all the data, but loads it in an easily accesible manner  
  %
  %   pic = PIC(h5FilePath)
  %   Bx = pic.Bx; % Bx is a (nt x nx x ny) matrix
  %   B = pic.B; % structure with 3 (nt x nx x ny) matrices  
  
  properties (SetAccess = immutable)
    software
    ndim
    file
    namelist
    particlebinning
    attributes
    info
    species
    parent
  end
  
  properties (Access = protected)
    % Access = protected – access from class or subclasses
    % Data can be arbitrary size, so the class contains a pointer to the 
    % data file and each time loads the data with
    fields_
    iteration_
    twpe_
    twci_
    xe_
    ye_
    ze_
    xi_
    yi_
    zi_
    grid_
    indices_
    it_
    ix_
    iy_
    iz_
%    wpewce_ = [];
%    mime_ = [];
    
  end
  
  properties (Dependent = true)
    % Can be checked when setting, for example, right size, right type
    % Dependent properties don't store a value and can't be assigned
    % a value in their set method.
    fields
    iteration
    twpe
    twci
    xe
    ye
    ze
    xi
    yi
    zi
    grid
    indices
    it
    ix
    iy
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
    charge
  end
      
  properties
    userData = []; % anything can be added here
    xevar
    zevar
    xivar
    zivar
  end
  
  methods
    function obj = PIC3D(h5filePath,nameList,particleBinningPath)
      % pic = PIC(pathFields) % michaels code
      % pic = PIC(pathFields,pathNameList,pathParticleBinning) - % smilei code
      %   pathParticleBinning should contain wildcard: ParticleBinning*.h5
      %
      % Fields:
      %   file - path to hdf5 file
      %   info - hdf5 structure and contents: h5info(file)
      %   mass - mass of species
      %   mime - mp/me
      %   teti - Te0/Tp0 - Harris sheet population
      %   wpewce - wpe0/wce0 (based on B0 and n0) 
      %   iteration - simulation iteration, one iteration is one step
      %   twpe0 - time in inverse electron plasma frequencies
      %   twci0 - time in inverse ion cyclotron frequencies
      %   xe, ze - x, z coordinate in terms of electron inertial lengths
      %   xi, zi - x, z coordinate in terms of proton inertial lengths
      %   m                   
      tic;
      obj.file = h5filePath; 
      obj.info = h5info(h5filePath);
      
      % Check if it's SMILEI or micPIC (mic stands for Michael)      
      obj.software = get_software(obj);
      
      if strcmp(obj.software,'Smilei') 
        if nargin == 1
          error('Wrong number of input arguments for Smilei simulation.')
        end
        if nargin >= 2
          % meta data
          obj.namelist = nameList;
          namelist = parse_namelist(obj);
          obj.species = namelist.name;
          obj.charge = namelist.charge;
          obj.mass = namelist.mass;
          obj.mime = namelist.mime;
          obj.wpewce = namelist.wpewce;
          obj.teti = namelist.teti;
          % grid
          obj.xe = namelist.xe;
          obj.ye = namelist.ye;
          obj.ze = namelist.ze;
          obj.attributes.particles_per_cell = namelist.particles_per_cell;          
        end
        if nargin == 3 % also load particle binning info
          obj.particlebinning = particleBinningPath;
          obj.attributes.deposited_quantity = namelist.deposited_quantity;
          obj.attributes.deposited_species = namelist.deposited_species;
        end
      elseif strcmp(obj.software,'micPIC')
        % meta data
        obj.charge = obj.get_charge;
        obj.mass = obj.get_mass;
        uniqueMass = sort(unique(obj.mass));
        obj.mime = uniqueMass(2)/uniqueMass(1); % second lightest/lightest
        obj.wpewce = h5read(h5filePath,'/simulation_information/wpewce');
        obj.teti = h5read(h5filePath,'/simulation_information/teti');
        % grid
        obj.xe = h5read(h5filePath,'/simulation_information/xe'); % de
        obj.ze = h5read(h5filePath,'/simulation_information/ze'); % de
        
        % The grid is offset for different variables (see figure 
        % micPIC_grid_layout in git repository). Needs to be used
        % when calculating derivatives (gradients, curls).

        % electorn intertial lengths
        dxe = (obj.xe(2) - obj.xe(1))/2;
        dze = (obj.ze(2) - obj.ze(1))/2;
        obj.xevar.Ey = obj.xe; % Ey on grid
        obj.zevar.Ey = obj.ze; % Ey on grid
        obj.xevar.Ex = obj.xe + dxe; % Ex offset half step to the right
        obj.zevar.Ex = obj.ze;       % Ex offset half step to the right
        obj.xevar.Bz = obj.xe + dxe; % Bz offset half step to the right
        obj.zevar.Bz = obj.ze;       % Bz offset half step to the right
        obj.xevar.Ez = obj.xe;       % Ez offset half step to the top
        obj.zevar.Ez = obj.ze + dze; % Ez offset half step to the top
        obj.xevar.Bx = obj.xe;       % Bx offset half step to the top
        obj.zevar.Bx = obj.ze + dze; % Bx offset half step to the top  
        obj.xevar.By = obj.xe + dze; % By offset half step to the right and top
        obj.zevar.By = obj.ze + dze; % By offset half step to the right and top

        % ion intertial lengths
        obj.xivar.Ey = obj.xevar.Ey/sqrt(obj.mime);
        obj.zivar.Ey = obj.zevar.Ey/sqrt(obj.mime);                    
        obj.xivar.Ex = obj.xevar.Ex/sqrt(obj.mime);
        obj.zivar.Ex = obj.zevar.Ex/sqrt(obj.mime);
        obj.xivar.Bz = obj.xevar.Bz/sqrt(obj.mime);
        obj.zivar.Bz = obj.zevar.Bz/sqrt(obj.mime);
        obj.xivar.Ez = obj.xevar.Ez/sqrt(obj.mime);
        obj.zivar.Ez = obj.zevar.Ez/sqrt(obj.mime);
        obj.xivar.Bx = obj.xevar.Bx/sqrt(obj.mime);
        obj.zivar.Bx = obj.zevar.Bx/sqrt(obj.mime);
        obj.xivar.By = obj.xevar.By/sqrt(obj.mime);
        obj.zivar.By = obj.zevar.By/sqrt(obj.mime);
      end
      
      obj.xi = obj.xe/sqrt(obj.mime);
      obj.yi = obj.ye/sqrt(obj.mime);
      obj.zi = obj.ze/sqrt(obj.mime);
      obj.grid = {1:1:numel(obj.xe),1:1:numel(obj.ye),1:1:numel(obj.ze)}; % originally, complete grid
      obj.ix = 1:1:numel(obj.xe);
      obj.iy = 1:1:numel(obj.ye);
      obj.iz = 1:1:numel(obj.ze);
      
        
      obj.iteration = get_iterations(obj);      
      obj.twpe = get_twpe(obj);      
      obj.twci = obj.twpe/(obj.wpewce*obj.mime);
      obj.indices_ = 1:numel(obj.iteration);
      obj.it = 1:1:numel(obj.twpe);
      
      obj.fields_ = get_fields(obj);
      if strcmp(obj.software,'Smilei') && nargin == 3
        obj.fields_ = cat(1,obj.fields_(:),unique(namelist.deposited_quantity(:)));
      end
      toc
      obj.parent = obj;
      obj.attributes = get_attributes(obj);
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
          obj.it_ = builtin('subsref',obj.it,s);
          if numel(idx(1).subs) == 3 % time and two spatial indices
            s = substruct(idx(1).type,idx(1).subs(2));
            newgrid{1} = builtin('subsref',obj.grid{1},s); 
            obj.xe_ = builtin('subsref',obj.xe,s); 
            obj.xi_ = builtin('subsref',obj.xi,s); 
            obj.ix_ = builtin('subsref',obj.ix,s); 
            s = substruct(idx(1).type,idx(1).subs(3));
            newgrid{2} = builtin('subsref',obj.grid{2},s);
            obj.ze_ = builtin('subsref',obj.ze,s); 
            obj.zi_ = builtin('subsref',obj.zi,s); 
            obj.iz_ = builtin('subsref',obj.iz,s); 
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
    function obj = xgrid(obj,inds)
      % pic.XGRID Get subset of x indices.
      %   pic.XGRID(100) - 
      %   pic.XGRID([100 110]) - 
      %   pic.XGRID(10:1:110) - 
      %
      % See also: PIC.ZLIM, PIC.TWPELIM, PIC.TWCILIM      
      obj = obj.subset('x',inds);      
    end
    function obj = ygrid(obj,inds)
      % pic.ZGRID Get subset of z indices.
      %   pic.ZGRID(100) - 
      %   pic.ZGRID([100 110]) - 
      %   pic.ZGRID(10:1:110) - 
      %
      % See also: PIC.XLIM, PIC.TWPELIM, PIC.TWCILIM      
      obj = obj.subset('y',inds);
    end
    function obj = zgrid(obj,inds)
      % pic.ZGRID Get subset of z indices.
      %   pic.ZGRID(100) - 
      %   pic.ZGRID([100 110]) - 
      %   pic.ZGRID(10:1:110) - 
      %
      % See also: PIC.XLIM, PIC.TWPELIM, PIC.TWCILIM      
      obj = obj.subset('z',inds);
    end
    function obj = tgrid(obj,inds)
      % pic.ZGRID Get subset of z indices.
      %   pic.ZGRID(100) - 
      %   pic.ZGRID([100 110]) - 
      %   pic.ZGRID(10:1:110) - 
      %
      % See also: PIC.XLIM, PIC.TWPELIM, PIC.TWCILIM      
      obj = obj.subset('t',inds);
    end
    function obj = xlim(obj,value,varargin)
      % pic.XLIM Get subset of xi (x/di).
      %   pic.XLIM(100) - gives closest index
      %   pic.XLIM([100 110]) - gives all indiced within range
      %   pic.XLIM(10:1:110) - gives all indices that matches exactly 
      %
      % See also: PIC.ZLIM, PIC.TWPELIM, PIC.TWCILIM
      inds = obj.ind_from_lim(obj.xi_,value,varargin{:});
      obj = obj.subset('x',inds);      
    end
    function obj = ylim(obj,value,varargin)
      % pic.ZLIM Get subset of zi (z/di).
      %   pic.ZLIM(100) - gives closest index
      %   pic.ZLIM([100 110]) - gives all indiced within range
      %   pic.ZLIM(10:1:110) - gives all indices that matches exactly 
      %
      % See also: PIC.XLIM, PIC.TWPELIM, PIC.TWCILIM 
      inds = obj.ind_from_lim(obj.zi_,value,varargin{:});
      obj = obj.subset('y',inds);
    end
    function obj = zlim(obj,value,varargin)
      % pic.ZLIM Get subset of zi (z/di).
      %   pic.ZLIM(100) - gives closest index
      %   pic.ZLIM([100 110]) - gives all indiced within range
      %   pic.ZLIM(10:1:110) - gives all indices that matches exactly 
      %
      % See also: PIC.XLIM, PIC.TWPELIM, PIC.TWCILIM 
      inds = obj.ind_from_lim(obj.zi_,value,varargin{:});
      obj = obj.subset('z',inds);
    end
    function obj = tlim(obj,value,varargin)
      % pic.TLIM Get subset of t (twci). Same as pic.TWCILIM.
      %   pic.TLIM(100) - gives closest index
      %   pic.TLIM([100 110]) - gives all indiced within range
      %   pic.TLIM(10:1:110) - gives all indices that matches exactly 
      %
      % See also: PIC.XLIM, PIC.ZLIM, PIC.TWCILIM, PIC.TWPELIM 
      inds = obj.ind_from_lim(obj.twci,value,varargin{:});    
      obj = obj.subset('t',inds);     
    end
    function obj = twpelim(obj,value,varargin)
      % pic.TWPELIM Get subset of t (twpe).
      %   pic.TWPELIM(1000) - gives closest index
      %   pic.TWPELIM([1000 10000]) - gives all indiced within range
      %   pic.TWPELIM(1000:200:10000) - gives all indices that matches exactly 
      %
      % See also: PIC.XLIM, PIC.ZLIM, PIC.TWCILIM 
      inds = obj.ind_from_lim(obj.twpe_,value,varargin{:});    
      obj = obj.subset('t',inds);
%       obj.twpe_ = obj.twpe_(inds);
%       obj.twci_ = obj.twci_(inds);
%       obj.it_ = obj.it_(inds);
%       obj.indices_ = obj.indices_(inds);
%       obj.iteration_ = obj.iteration_(inds);
    end
    function obj = twcilim(obj,value,varargin)
      % pic.TWCILIM Get subset of t (twci).
      %   pic.TWCILIM(100) - gives closest index
      %   pic.TWCILIM([100 110]) - gives all indiced within range
      %   pic.TWCILIM(10:1:110) - gives all indices that matches exactly 
      %
      % See also: PIC.XLIM, PIC.ZLIM, PIC.TWPELIM 
      inds = obj.ind_from_lim(obj.twci,value,varargin{:});    
      obj = obj.subset('t',inds);
    end
    function obj = subset(obj,comp,inds)
      % select subset of indices
      switch comp
        case 't'
          obj.twpe_ = obj.twpe_(inds);
          obj.twci_ = obj.twci_(inds);
          obj.it_ = obj.it_(inds);
          obj.indices_ = obj.indices_(inds);
          obj.iteration_ = obj.iteration_(inds);
        case 'x'
          obj.xe_ = obj.xe_(inds);
          obj.xi_ = obj.xi_(inds);      
          obj.grid_{1} = obj.grid_{1}(inds);
          obj.ix_ = obj.grid_{1};
          if strcmp(obj.software,'Smilei') 
          
          else
            field_names = fields(obj.xivar);
            for ifield = 1:numel(field_names)
              obj.xivar.(field_names{ifield}) = obj.xivar.(field_names{ifield})(inds);
            end
          end
        case 'y'
          obj.ye_ = obj.ye_(inds);
          obj.yi_ = obj.yi_(inds);      
          obj.grid_{2} = obj.grid_{2}(inds);
          obj.iy_ = obj.grid_{2};
          if strcmp(obj.software,'Smilei') 
          
          else
            field_names = fields(obj.xivar);
            for ifield = 1:numel(field_names)
              obj.xivar.(field_names{ifield}) = obj.xivar.(field_names{ifield})(inds);
            end
          end
        case 'z'
          obj.ze_ = obj.ze_(inds);
          obj.zi_ = obj.zi_(inds);      
          obj.grid_{3} = obj.grid_{3}(inds);
          obj.iz_ = obj.grid_{3};  
          if strcmp(obj.software,'Smilei') 
          
          else
            field_names = fields(obj.zivar);
            for ifield = 1:numel(field_names)
              obj.zivar.(field_names{ifield}) = obj.zivar.(field_names{ifield})(inds);
            end
          end
      end
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
      doAdjustCMap = 0;
      
      nfields = numel(fields);
      ntimes = obj.length;
      
      have_options = 0;
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
          %contourf(hca,obj.xi,obj.zi,squeeze(data)',0:0.5:25)
          imagesc(hca,obj.xi,obj.zi,squeeze(data)')
          hb = colorbar('peer',hca);
          hb.YLabel.String = fields{ifield};
          if doAdjustCLim
            hca.CLim = clims{ifield};            
            %colormap(cmap)
          end
          if doAdjustCMap
            if isa(cmaps,'cell')
              hca.CLim = cmaps{ifield};
            elseif isnumeric(cmaps)
              colormap(hca,cmaps)
            end
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
        h(1).Title.String = ['t\omega_{pe} = ' sprintf('%g',obj.twpe(itime)) ', t\omega_{ci} = ' sprintf('%g',obj.twci(itime))];
        drawnow
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
    function [all_im, map] = movie(obj,varstrs,varargin)
      % pic.MOVIE(obj,varstrs,varargin) Makes movie and gif.
      %   pic.MOVIE(varstrs,'inp1',arg1,'inp2',arg2,...)
      %   
      %   Miinum input:
      %     varstrs - NxM cell array with variables or expressions with variables (see pic.get_exp)
      %   Additional input:
      %     'A' - plots A contours with given interval/step in A
      %     'cmap' - NxM cell array of colormaps
      %     'clim' - NxM cell array of clims
      %     'filename' - filename and path (excluding suffix/extension)
      %     'cbarlabels' - NxM cell array with colorbar labels (default is varstrs text)
      %     'traj' - PICTraj object, adds trajectories of object in each frame
      %     'trajcolordot' - color of trajectories/dots, can be:
      %       1x3 numeric array
      %       Ntrajx3 numeric array      
      %       string with PICTraj quantity as color (only implemented for scalars I think, e.g. 'x0')
      %
      % Examples:
      %   pic.MOVIE({'Ey'})
      %   pic.zlim([0 25]*0.99).MOVIE({'tzz(4)'},'A',1,'cmap',pic_colors('waterfall'),'clim',{[0 0.015]},'filename',filename);
      %
      % See also: PIC.GET_EXP
      
      % Subsref error: Does not accept two outputs
      
      % Default options, values
      doCBarLabels = 0;
      doVideo = 1;
      doGif = 1;
      doGifBackLoop = 1;
      doAdjustCLim = 0;
      cmap = pic_colors('blue_red');
      doA = 0;
      colorA = [0 0 0];
      doAdjustCMap = 0;
      fileName = 'movie'; % .mp4/.gif added later
      doTrajectories = 0;
      colorTrajDot = [0 0 0];
      doColorTrajDot = 0;
      doTrajTail = 1;
      ntTail = 4;
      stepA = 0.5;
      trajargs = {};
      print_flag = 1;
      doSmooth = 0;
      npSmooth = 1;
      colorMapTraj = pic_colors('waterfall');
      
      ntimes = obj.length;
      
      have_options = 0;
      nargs = numel(varargin);      
      if nargs > 0, have_options = 1; args = varargin(:); end
      
      while have_options
        l = 1;
        %lower(args{1})
        switch(lower(args{1}))
          case {'traj','trajectories','tr','orbits'}
            doTrajectories = 1;
            tr = args{2};
            l = 2;         
          case {'trajcolordot'}
            doColorTrajDot = 1;
            cTrajDot = args{2};
            l = 2;
          case 'a'
            doA = 1;
            stepA = args{2};
            l = 2;
          case {'cola'}
            colorA = args{2};
            l = 2;
          case 'smooth'
            doSmooth = 1;
            npSmooth = args{2};
            l = 2;
          case {'trajargs','trarg','argtraj','argtr'}
            trajargs = args{2};
            l = 2;
          case 'cbarlabels'
            doCBarLabels = 1;
            cBarLabels = args{2};
            l = 2;
          case 'clim'
            l = 2;
            doAdjustCLim = 1;  
            clims = args{2};
          case 'cmap'
            l = 2;
            doAdjustCMap = 1;
            cmaps = args{2};
          case {'path','filename'}
            l = 2;
            fileName = args{2};
          case 'gif'
            l = 1;
            doGif = 1;
            doVideo = 0;            
          otherwise
            warning(sprintf('Input ''%s'' not recognized.',args{1}))            
        end        
        args = args(l+1:end);
        if isempty(args), break, end    
      end
         
      if exist([fileName '.mp4'],'file')
        question = 'File already exists, do you want to overwrite it? [1/0] >';
        print_flag = irf_ask(question,[],1); 
      end
      if print_flag == 0
        disp('Not making movie.')
        return;
      end

      if doColorTrajDot % is this needed here?
        if ischar(cTrajDot)
          colorTrajDot = cTrajDot; % color will be assigned further down, when trajectory is loaded
        elseif 0%any(or(size(cTrajDot)==[tr.ntr 1],size(cTrajDot)==[1 tr.ntr]))
          if numel(unique(cTrajDot))<=7
             % matlab colors
            colors = [     0    0.4470    0.7410;
                      0.8500    0.3250    0.0980;
                      0.9290    0.6940    0.1250;
                      0.4940    0.1840    0.5560;
                      0.4660    0.6740    0.1880;
                      0.3010    0.7450    0.9330;
                      0.6350    0.0780    0.1840];
            for itr = 1:tr.ntr
              colorTrajDot(itr,:) = colors(fix(cTrajDot(itr)),:);
            end
          else
            colorTrajDot = reshape(cTrajDot,tr.ntr,1);
          end
        else
          colorTrajDot = cTrajDot
        end
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
          h(irow,icol).Position(2) = h(irow,icol).Position(2)+0.05;
          h(irow,icol).Position(4) = h(irow,icol).Position(4)-0.05/nrows;
        end
      end
      
      if doVideo
        vidObj = VideoWriter([fileName '.mp4'],'MPEG-4');
        vidObj.FrameRate = 10;
        open(vidObj);        
      end
      if doGif
        iframe = 0;
      end
      
      disp('Adjust figure size, then hit any key to continue.')
      pause
      
      for itime = 1:obj.nt
        tmp_obj = obj.twcilim(obj.twci(itime));
        hleg = gobjects(0);
        ivar = 0;
        for irow = 1:nrows
          for icol = 1:ncols
            ivar = ivar + 1;
            isHoldOn = 0;
            ip = sub2ind([ncols nrows],icol,irow);
            hca = h(ip);

            % check if input demand som andditional input, e.g. n(1)
            if strfind(varstrs{ivar},'Acont') % only plot A contours
              plot(hca,NaN,NaN)
              doA = 1;
            else              
%               if strfind(varstrs{ivar},'(')
%                 ind1 = strfind(varstrs{ivar},'(');
%                 ind2 = strfind(varstrs{ivar},')');
%                 %varsplit = regexp(varstrs{ivar}, '(?<var>\w+)\W+(?<ind>\d)\W+','names');
%                 indstr = varstrs{ivar}(ind1+1:ind2-1);
%                 varstr =  varstrs{ivar}(1:ind1-1);
%                 var = tmp_obj.(varstr)(eval(indstr));
%               else
%                 var = tmp_obj.(varstrs{ivar});
%               end        
              var = tmp_obj.get_exp(varstrs{ivar});              
              if doSmooth
                var = smooth2(var,npSmooth);
              end
              imagesc(hca,tmp_obj.xi,tmp_obj.zi,var');
              %contourf(hca,tmp_obj.xi,tmp_obj.zi,var',0:0.25:25);
              %contourf(hca,tmp_obj.xi,tmp_obj.zi,var',2:0.2:8);
              hb(ivar) = colorbar('peer',hca);
              if doCBarLabels
                hb(ivar).YLabel.String = cBarLabels{ivar};
              else
                hb(ivar).YLabel.String = varstrs{ivar};                
              end
              if doAdjustCLim
                hca.CLim = clims{ip};
                %colormap(cmap)
              end
              if doAdjustCMap
              if isa(cmaps,'cell')
                colormap(hca,cmaps{ip});
              elseif isnumeric(cmaps)
                colormap(hca,cmaps)
              end
              %colormap(cmap)
            end
            end
            hca.XLabel.String = 'x (d_i)';
            hca.YLabel.String = 'z (d_i)';
            hca.YDir = 'normal';
            if doA
              clim = hca.CLim;
              A = tmp_obj.A;
              %stepA = 1;
              levA = floor(min(A(:))/stepA)*stepA:stepA:ceil(max(A(:))/stepA)*stepA;
              %levA = stepA;
              iAx = 1:5:obj.nx;
              iAz = 1:5:obj.nz;
              hold(hca,'on')
              contour(hca,tmp_obj.xi(iAx),tmp_obj.zi(iAz),A(iAx,iAz)',levA,'color',colorA)
              hold(hca,'off')
              hca.CLim = clim; 
            end        
            if doTrajectories
              hold(hca,'on')
              for itr = 1:tr.ntr
                tt = tr(itr).t;
                xx = tr(itr).x;
                zz = tr(itr).z;
                idup = find(diff(tr(itr).t)==0);
                tt(idup) = [];
                xx(idup) = [];
                zz(idup) = [];
                if tmp_obj.twci == tt(1)
                  xnow = xx(1);
                  znow = zz(1); 
                  xtail = NaN;
                  ztail = NaN;
                elseif tmp_obj.twci == tt(end)
                  xnow = xx(end);
                  znow = zz(end);
                  xtail = NaN;
                  ztail = NaN;
                else
                  if 1
                    xnow = interp1(tt,xx,tmp_obj.twci);
                    znow = interp1(tt,zz,tmp_obj.twci);
%                   else
%                     ind_ = itime + [-5:0];
%                     ind_ = ind_(ind_ > 0);
%                     tmp_obj_ = obj.twcilim(obj.twci(ind_));
%                     xnow = interp1(tt,xx,tmp_obj_.twci);
%                     znow = interp1(tt,zz,tmp_obj_.twci);
                  end
                  if doTrajTail
                    if itime == 1
                      xtail = NaN;
                      ztail = NaN;
                    else
                      ind_ = itime + [-ntTail:0];
                      ind_ = ind_(ind_ > 0);
                      it1_tail = find(tt>obj.twci(ind_(1)),1,'first');
                      it2_tail = find(tt<obj.twci(ind_(end)),1,'last');
                      it_tail = it1_tail:it2_tail;
                      xtail = xx(it_tail);                    
                      ztail = zz(it_tail);
                    end
                  end
                end
                %plot(hca,tr(itr).x,tr(itr).z,'k')
                if doColorTrajDot
                  if ischar(colorTrajDot)
                    var = tr(itr).(colorTrajDot); 
                    var(idup) = [];
                    coldot = interp1(tt,var,tmp_obj.twci);
                  else
                    if numel(colorTrajDot) == 3
                      coldot = colorTrajDot;
                    elseif numel(colorTrajDot) == tr.ntr
                      %coldot = colorTrajDot(itr,:);
                      crange = (cTrajDot(itr)-0.99999*min(cTrajDot))/(max(cTrajDot)-0.9999*min(cTrajDot));
                      coldot = interp1(0:size(colorMapTraj,1)-1,colorMapTraj,crange*(size(colorMapTraj,1)-1));
                    elseif size(colorTrajDot) == [tr.ntr,3]
                      coldot = colorTrajDot(itr,:);
                    end
                  end
                else
                  coldot = [0 0 0];
                end
                %plot(hca,xnow,znow,'color',coldot,'markerSize',20,'marker','.','linestyle','none')
                sc = scatter(hca,xnow,znow,30,coldot,'Marker','o','MarkerFaceColor','flat','MarkerEdgeColor',[0 0 0],trajargs{:});
                %hp = plot(hca,tr(itr).x0,tr(itr).z0,'ko');
                if doTrajTail
                 plot(hca,xtail,ztail,'color',coldot)  
                end
              end
              hold(hca,'off')
            end
            hca.XLim = obj.xi([1 end]);
            hca.YLim = obj.zi([1 end]);
            hca.FontSize = 14;
           % drawnow;
          end
        end 
        %drawnow;
        h(1).Title.String = sprintf('twpe = %.0f, twci = %.1f',tmp_obj.twpe,tmp_obj.twci);
        compact_panels(0.01)
%        h(3).Position(3) = h(2).Position(3);
        % Collect frames
        pause(1)
        if doVideo
          set(gcf,'color','white');
          currFrame = getframe(gcf);
          writeVideo(vidObj,currFrame);
        end
        if doGif
          if 1 % collect frames, for making gif
            iframe = iframe + 1;    
            nframes = obj.nt;
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
      end
      
      % Write gif
      if doGif
        imwrite(all_im,map,[fileName,'.gif'],'DelayTime',0,'LoopCount',inf)
      end
      if doGif && doGifBackLoop
        imwrite(cat(4,all_im,all_im(:,:,:,end:-1:1)),map,[fileName,'_loopback.gif'],'DelayTime',0,'LoopCount',inf)              
      end
      
      %hlinks = linkprop(h,{'XLim','YLim'});
      %set(gcf,'userdata',{'hlinks',hlinks})
      if nargout == 1
        varargout{1} = h;
      elseif nargout == 2
        varargout{1} = h;
        varargout{2} = hb;
      elseif nargout == 3
        varargout{1} = h;
        varargout{2} = hb;
        varargout{3} = hlinks;
      end
                     
    end
    function varargout = movie_line(obj,dim,varstrs_all,varargin)
      % Plots variables directly loaded from file
      % h = pic.plotline(dim,{{'Ex','Ey','Ez'},{'Bx','By','Bx'}})
      % Plot variables as a function of dim: 'x' or 'z'
      
      doVideo = 1;
      doGif = 0;
      doYLim = 0;
      fileName = 'movie'; % deafult filenam (including path, mp4 added later)
      doSmooth = 0;
      doColor = 0;
      doVarLabels = 0;
     
      % Which dimension to plot against
      if strcmp(dim,'x')
        plot_dep = obj.xi;
        dep_lim = obj.xi([1 end]);
        sum_dim = 2;
        range = obj.zi([1 end]);
        rangestr = 'z';
      elseif strcmp(dim,'z')
        plot_dep = obj.zi;
        dep_lim = obj.zi([1 end]);
        sum_dim = 1;
        range = obj.xi([1 end]);
        rangestr = 'x';        
      else
        error(sprintf('Unknown dependent dimension %s. Must be ''x'' or ''z''.',dim))
      end
            
      % Check additional input      
      have_options = 0;
      nargs = numel(varargin);      
      if nargs > 0, have_options = 1; args = varargin(:); end      
      while have_options
        l = 1;
        switch(lower(args{1}))
          case 'ylim'
            l = 2;
            doYLim = 1;  
            ylims = args{2};
          case 'color'
            l = 2;
            doYLim = 1;  
            ylims = args{2};
          case {'path','filename'}
            l = 2;
            fileName = args{2};
          case {'smooth'}
            l = 2;
            doSmooth = 1;
            npSmooth = args{2};
          otherwise 
            warning(sprintf('Unknown argument %s.',args{1}))
        end
        args = args(l+1:end);  
        if isempty(args), break, end    
      end
      
      [nrows,ncols] = size(varstrs_all);    
      npanels = nrows*ncols;
      ip = 0;
      for irow = 1:nrows
        for icol = 1:ncols
          ip = ip + 1;
          h(irow,icol) = subplot(nrows,ncols,ip);
        end
      end
      
        
      if doVideo
        vidObj = VideoWriter([fileName '.mp4'],'MPEG-4');
        open(vidObj);
      end
      if doGif
        iframe = 0;
      end
      
      disp('Adjust figure size, then hit any key to continue.')
      pause
      
      for it = 1:obj.nt
        obj_tmp = obj.subset('t',it);
        hleg = gobjects(0);
        ip = 0;
        for irow = 1:nrows
          for icol = 1:ncols
            ip = ip + 1;
            isHoldOn = 0;
            %ip = sub2ind([ncols nrows],icol,irow);
            hca = h(irow,icol);
            % check if input demand som andditional input, e.g. n(1)
            varstrs = varstrs_all{ip};
            nvars = numel(varstrs);
            for ivar = 1:nvars
              varstr = varstrs{ivar};
              if strcmp(varstr(1),'-')
                varmult = -1;
                varstr = varstr(2:end);
              else
                varmult = 1;
              end
              if strfind(varstr,'(')
                ind1 = strfind(varstr,'(');
                ind2 = strfind(varstr,')');            
                indstr = varstr(ind1+1:ind2-1);
                varstr =  varstr(1:ind1-1);
                var = varmult*obj_tmp.(varstr)(eval(indstr));
              else
                var = varmult*obj_tmp.(varstr);
              end
              displayname = [varstrs{ivar}];
              
              var = squeeze(mean(var,sum_dim));
              if doSmooth
                var = smooth(var,npSmooth);
              end
              plot(hca,plot_dep,var,'DisplayName',displayname);
              if doYLim
                hca.YLim = ylims{ip};
              end
              if not(isHoldOn)
                hold(hca,'on')
                isHoldOn = 1;
              end              
              hca.XLabel.String = [dim ' (d_i)'];
            end
            hold(hca,'off')
            hca.XGrid = 'on';
            hca.YGrid = 'on';
            hleg(ip) = legend(hca,'location','eastoutside'); % this may make the panels of different widths, fix below
          end
        end
        drawnow;
        leftpos = 1;
        for ip = 1:npanels, panel_width(ip) = h(ip).Position(3); end      
        for ip = 1:npanels, h(ip).Position(3) = min(panel_width); end        
        titlestring = sprintf('twpe = %g, twci = %.1f, %s = [%.2f %.2f]',obj_tmp.twpe,obj_tmp.twci,rangestr,range(1),range(2));
        h(1).Title.String = titlestring;
        compact_panels(0.01)
        hlinks = linkprop(h,{'XLim'});
        h(1).XLim = dep_lim;
        set(gcf,'userdata',{'hlinks',hlinks})
            
        % Collect frames
        pause(1)
        if doVideo
          set(gcf,'color','white');
          currFrame = getframe(gcf);
          writeVideo(vidObj,currFrame);
        end
        if doGif
          if 1 % collect frames, for making gif
            iframe = iframe + 1;    
            nframes = pic0.nt;
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
      end
      
      if nargout == 1
        varargout{1} = h;
      elseif nargout == 2
        varargout{1} = h;
        varargout{2} = hleg;
      elseif nargout == 2
        varargout{1} = h;
        varargout{2} = hleg;
        varargout{3} = hlinks;
      end
    end
    function varargout = plot_map(obj,varargin)
      % pic.PLOTMAP Plots variables in (x,z) map
      %   h = pic.PLOTMAP(varstrs,'inp1','arg1',...);
      %   [h,hb] = pic.PLOTMAP(varstrs,'inp1','arg1',...);
      %   [h,hb,hlinks] = pic.PLOTMAP(varstrs,'inp1','arg1',...);
      %
      % Output:
      %   h - handle to axes
      %   hb - handle to colorbars
      %   hlinks - handle to linked stuff, see linkprop.
      
     
      
      % Defaults
      doCBarLabels = 0;
      doSep = 0;
      doA = 0;
      doAdjustCLim = 0;
      cmap = pic_colors('blue_red');      
      doAdjustCMap = 0;
      doSmooth = 0;
      plotInAxes = 0;
      doCompact = 1;
      
      % Check input
      % Check for axes
      [ax,args,nargs] = irf.axescheck(varargin{:}); 
      if not(isempty(ax)); plotInAxes = 1; doCompact = 0; end
      
      varstrs = args{1};
      args = args(2:end);
      nargs = nargs-1;
      have_options = 0;      
      if nargs > 0, have_options = 1; end
      
      while have_options
        l = 1;
        switch(lower(args{1}))
          case 'a'
            doA = 1;
            stepA = args{2};
            l = 2;
          case 'clim'
            l = 2;
            doAdjustCLim = 1;  
            clims = args{2};
          case 'cmap'
            l = 2;
            doAdjustCMap = 1;
            cmaps = args{2};
          case 'smooth'
            l = 2;
            doSmooth = 1;
            npSmooth = args{2};            
          case 'sep'
            doSep = 1;
            l = 1;
          case 'cbarlabels'
            doCBarLabels = 1;
            cBarLabels = args{2};
            l = 2;
          otherwise 
            warning(sprintf('Unknown argument %s.',args{1}))
        end
        args = args(l+1:end);  
        if isempty(args), break, end    
      end
      
      if doA
        A = obj.A;
        if numel(stepA) == 1
          levA = floor(min(A(:))/stepA)*stepA:stepA:ceil(max(A(:))/stepA)*stepA;
        else
          levA = stepA;
        end
        iAx = 1:5:obj.nx;
        iAz = 1:5:obj.nz;
      end
      if doSep
        if not(doA)
          A = obj.A;
        end        
        [saddle_locations,saddle_values] = saddle(A,'sort');
        sepA = saddle_values(1);
      end
      
      [nrows,ncols] = size(varstrs);      
      npanels = nrows*ncols;
      if plotInAxes && numel(ax) == npanels
        h = ax;
      else
        ip = 0;
        for irow = 1:nrows
          for icol = 1:ncols
            ip = ip + 1;
            h(irow,icol) = subplot(nrows,ncols,ip);
          end
        end
      end
      
      % Plot
      hleg = gobjects(0);
      ivar = 0;
      for irow = 1:nrows
        for icol = 1:ncols
          ivar = ivar + 1;
          isHoldOn = 0;
          ip = sub2ind([ncols nrows],icol,irow);
          hca = h(ip);
        
          % check if input demand som andditional input, e.g. n(1)

%           if strfind(varstrs{ivar},'(')
%             ind1 = strfind(varstrs{ivar},'(');
%             ind2 = strfind(varstrs{ivar},')');
%             %varsplit = regexp(varstrs{ivar}, '(?<var>\w+)\W+(?<ind>\d)\W+','names');
%             indstr = varstrs{ivar}(ind1+1:ind2-1);
%             varstr =  varstrs{ivar}(1:ind1-1);
%             var = obj.(varstr)(eval(indstr));
%           else
%             var = obj.(varstrs{ivar});
%           end  
          var = obj.get_exp(varstrs{ivar}); 
          if doSmooth
            var = smooth2(var,npSmooth);
          end
          imagesc(hca,obj.xi,obj.zi,var');
          hb(ivar) = colorbar('peer',hca);
          hb(ivar).YLabel.String = varstrs{ivar};
          hca.XLabel.String = 'x (d_i)';
          hca.YLabel.String = 'z (d_i)';
          hca.YDir = 'normal';
          clim = hca.CLim;
          if doA
            hold(hca,'on')
            contour(hca,obj.xi(iAx),obj.zi(iAz),A(iAx,iAz)',levA,'k')
            hold(hca,'off')
          end
          if doSep
            hold(hca,'on')
            contour(hca,obj.xi(iAx),obj.zi(iAz),A(iAx,iAz)',sepA*[1 1],'k--')
            hold(hca,'off')
          end
          hca.CLim = clim;
          if doAdjustCLim
            hca.CLim = clims{ip};
            %colormap(cmap)
          end
          if doAdjustCMap
          if isa(cmaps,'cell')
            colormap(hca,cmaps{ip});
          elseif isnumeric(cmaps)
            colormap(hca,cmaps)
          end
          %colormap(cmap)
        end
          if doCBarLabels
            hb(ivar).YLabel.String = cBarLabels{ivar};
          else
            hb(ivar).YLabel.String = varstrs{ivar};                
          end
          %hb(ivar).FontSize = 14;
         % drawnow;
        end
      end 
      %drawnow;
      h(1).Title.String = sprintf('twpe = %g, twci = %g',obj.twpe,obj.twci);
      %doCompact
      if doCompact, compact_panels(0.01); end
      hlinks = linkprop(h(:),{'XLim','YLim'});
      set(gcf,'userdata',{'hlinks',hlinks})
      if nargout == 1
        varargout{1} = h;
      elseif nargout == 2
        varargout{1} = h;
        varargout{2} = hb;
      elseif nargout == 3
        varargout{1} = h;
        varargout{2} = hb;
        varargout{3} = hlinks;
      end
    end
    function varargout = plot_timemap(obj,dim,varstrs,varargin)
      % Plots variables directly loaded from file
      % h = pic.PLOTTIMEMAP('zt',{'Ey','Ez'})
      % h = pic.PLOTTIMEMAP('tx',{'Ey','Ez'},'A')
      
      % Defaults
      doA = 0;
      doXline = 0;
      
      % Which dimension to plot against
      if strfind(dim,'x')
        if strfind(dim,'x') == 1 % x on x-axis
          plot_depx = obj.xi;
          plot_depy = obj.twci;
          permuteorder = [2 1];
        else
          plot_depx = obj.twci;
          plot_depy = obj.xi;
          permuteorder = [1 2];
        end                
        sum_dim = 2;
        range = obj.zi([1 end]);
        rangestr = 'z';
      elseif strfind(dim,'z')
        if strfind(dim,'z') == 1 % x on x-axis
          plot_depx = obj.zi;
          plot_depy = obj.twci;
          permuteorder = [2 1];
        else
          plot_depx = obj.twci;
          plot_depy = obj.zi;
          permuteorder = [1 2];
        end                
        sum_dim = 1;
        range = obj.xi([1 end]);
        rangestr = 'x';        
      else
        error(sprintf('Unknown dependent dimension %s. Must be ''xt'', ''tx'', ''zt'', or ''tz''.',dim))
      end
      lim_depx = plot_depx([1 end]);
      lim_depy = plot_depy([1 end]);
      
%       [ax,args,nargs] = irf.axescheck(varargin{:}); 
%       if not(isempty(ax)); plotInAxes = 1; doCompact = 0; end
%       
%       varstrs = args{1};
%       args = args(2:end);
%       nargs = nargs-1;
%       have_options = 0;      
%       if nargs > 0, have_options = 1; end
      
      % Check input
      have_options = 0;
      nargs = numel(varargin);      
      if nargs > 0, have_options = 1; args = varargin(:); end      
      while have_options
        l = 1;
        switch(lower(args{1}))
          case 'a'
            doA = 1;
            stepA = args{2};
            l = 2;
          case 'xline'
            doXline = 1;
            l = 1;
          otherwise 
            warning(sprintf('Unknown argument %s.',args{1}))
        end
        args = args(l+1:end);  
        if isempty(args), break, end    
      end
      
      if doA
        A = obj.A;
        A = squeeze(mean(A,sum_dim));
        A = permute(A,permuteorder);
        %stepA = 1;
        levA = floor(min(A(:))/stepA)*stepA:stepA:ceil(max(A(:))/stepA)*stepA;
        if strfind(dim,'t') == 1
          iAdepx = 1:1:numel(plot_depx);
          iAdepy = 1:5:numel(plot_depy);
        else
          iAdepx = 1:5:numel(plot_depx);
          iAdepy = 1:1:numel(plot_depy);
        end
      end
      if doXline
        switch dim
          case 'xt'
            x_xline = pic.z_xline;
            y_xline = pic.twci;
          case 'zt'
            x_xline = pic.z_xline;
            y_xline = pic.twci;
          case 'tx'
            y_xline = pic.z_xline;
            x_xline = pic.twci;
          case 'tz'
            y_xline = pic.z_xline;
            x_xline = pic.twci;
        end
      end
      [nrows,ncols] = size(varstrs);      
      for ip = 1:nrows*ncols
        h(ip) = subplot(nrows,ncols,ip);
      end
      hb = gobjects(0);
      for ivar = 1:nrows*ncols
        hca = h(ivar);
        % check if input demand som andditional input, e.g. n(1)
        tic
%         if strfind(varstrs{ivar},'(')
%           ind1 = strfind(varstrs{ivar},'(');
%           ind2 = strfind(varstrs{ivar},')');
%           %varsplit = regexp(varstrs{ivar}, '(?<var>\w+)\W+(?<ind>\d)\W+','names');
%           indstr = varstrs{ivar}(ind1+1:ind2-1);
%           varstr =  varstrs{ivar}(1:ind1-1);
%           var = obj.(varstr)(eval(indstr));
%         else
%           var = obj.(varstrs{ivar});
%         end
        var = obj.get_exp(varstrs{ivar});
        toc
        var = permute(squeeze(nanmean(var,sum_dim)),permuteorder);
        pcolor(hca,plot_depx,plot_depy,var);
        shading(hca,'flat')
        hb(ivar) = colorbar('peer',hca);
        hb(ivar).YLabel.String = varstrs{ivar};
        hca.XLabel.String = [dim(1) ' (d_i)'];
        hca.YLabel.String = [dim(2) ' (d_i)'];
        hca.YDir = 'normal';
        clim = hca.CLim;
        if doA
          hold(hca,'on')
          contour(hca,plot_depx(iAdepx),plot_depy(iAdepy),A(iAdepy,iAdepx),levA,'k')
          hold(hca,'off')
        end
        if doXline
          hold(hca,'on')
          plot(hca,x_xline,y_xline,'linewidth',1,'color',[0 0 0])
          hold(hca,'off')
        end
        hca.CLim = clim;   
        drawnow;
      end
      drawnow;
      h(1).Title.String = sprintf('%s = [%.2f %.2f]',rangestr,range(1),range(2));
      compact_panels(0.01)
      hlinks = linkprop(h,{'XLim','YLim'});
      set(gcf,'userdata',{'hlinks',hlinks})
      if nargout == 1
        varargout{1} = h;
      elseif nargout == 2
        varargout{1} = h;
        varargout{2} = hb;
      elseif nargout == 3
        varargout{1} = h;
        varargout{2} = hb;
        varargout{3} = hlinks;
      end
    end
    function varargout = plot_line(obj,varargin)
     % function varargout = plot_line(obj,dim,varstrs_all,varargin)
      % Plots variables directly loaded from file
      % h = pic.plotline(dim,{{'Ex','Ey','Ez'},{'Bx','By','Bx'}})
      % Plot variables as a function of dim: 'x' or 'z'
       
      doSmooth = 0;
      doVertical = 0;
      plotInAxes = 0;
      doCompact = 1;
                  
      % Check additional input     
      [ax,args,nargs] = irf.axescheck(varargin{:});       
      if not(isempty(ax)); plotInAxes = 1; doCompact = 0; end
      dim = args{1}; % required input
      varstrs_all = args{2}; % required input
      args = args(3:end);
      nargs = nargs-2;
      have_options = 0;      
            
      if nargs > 0, have_options = 1; end      
      while have_options
        l = 1;
        switch(lower(args{1}))
          case 'smooth'
            doSmooth = 1;
            npSmooth = args{2};
            l = 2;
          case {'vertical','vert'}
            doVertical = 1;
            l = 1;
          otherwise 
            warning(sprintf('Unknown argument %s.',args{1}))
        end
        args = args(l+1:end);  
        if isempty(args), break, end    
      end
      
      
      % Which dimension to plot against
      if strcmp(dim,'x')
        plot_dep = obj.xi;
        dep_lim = obj.xi([1 end]);
        sum_dim = 2;
        range = obj.zi([1 end]);
        rangestr = 'z';
      elseif strcmp(dim,'z')
        plot_dep = obj.zi;
        dep_lim = obj.zi([1 end]);
        sum_dim = 1;
        range = obj.xi([1 end]);
        rangestr = 'x';        
      else
        error(sprintf('Unknown dependent dimension %s. Must be ''x'' or ''z''.',dim))
      end
      
      [nrows,ncols] = size(varstrs_all);    
      npanels = nrows*ncols;
      if plotInAxes && numel(ax) == npanels
        h = ax;
      else
        ip = 0;
        for irow = 1:nrows
          for icol = 1:ncols
            ip = ip + 1;
            h(irow,icol) = subplot(nrows,ncols,ip);
          end
        end
      end
      
      hleg = gobjects(0);
      ip = 0;
      for irow = 1:nrows % Plotting
        for icol = 1:ncols
          ip = ip + 1;
          isHoldOn = 0;
          %ip = sub2ind([ncols nrows],icol,irow);
          hca = h(irow,icol);
          % check if input demand som andditional input, e.g. n(1)
          varstrs = varstrs_all{ip};
          nvars = numel(varstrs);
          for ivar = 1:nvars
            varstr = varstrs{ivar};
            var = obj.get_exp(varstr);
            if ivar == 2
              hold(hca,'on')
            end
            var = squeeze(mean(var,sum_dim));
            ntimes = obj.nt;size(var,2);
            for itime = 1:ntimes
              displayname = [varstrs{ivar} ' (t\omega_{pe}=' sprintf('%.0f)',obj.twpe(itime))];
              if ntimes == 1
                plot_var = var;                              
              else
                plot_var = var(:,itime);
              end
              if doSmooth 
                plot_var = smooth(plot_var,npSmooth);
              end
              if doVertical
                plot(hca,plot_var,plot_dep,'DisplayName',displayname);
                hca.YLabel.String = [dim ' (d_i)'];
              else
                plot(hca,plot_dep,plot_var,'DisplayName',displayname);
                hca.XLabel.String = [dim ' (d_i)'];
              end
              if not(isHoldOn)
                hold(hca,'on')
                isHoldOn = 1;
              end
            end
            
          end
          hold(hca,'off')
          hca.XGrid = 'on';
          hca.YGrid = 'on';
          if doVertical
            hleg(ip) = legend(hca,'location','northoutside'); % this may make the panels of different widths, fix below
          else
            hleg(ip) = legend(hca,'location','eastoutside'); % this may make the panels of different widths, fix below
          end
        end
      end
      drawnow;      
      leftpos = 1;
      if doVertical
        for ip = 1:npanels, panel_width(ip) = h(ip).Position(4); end      
        for ip = 1:npanels, h(ip).Position(4) = min(panel_width); end        
      else
        for ip = 1:npanels, panel_width(ip) = h(ip).Position(3); end      
        for ip = 1:npanels, h(ip).Position(3) = min(panel_width); end
      end
      if ntimes == 1
        titlestring = sprintf('twpe = %g, twci = %g, %s = [%.2f %.2f]',obj.twpe,obj.twci,rangestr,range(1),range(2));
      else
        titlestring = sprintf('twpe = %g-%g, twci = %g-%g, %s = [%.2f %.2f]',obj.twpe(1),obj.twpe(end),obj.twci(1),obj.twci(end),rangestr,range(1),range(2));
      end
      if doSmooth
        titlestring = [titlestring ', n_p^{smooth} = ' num2str(npSmooth)];
      end
      h(1).Title.String = titlestring;
      if doCompact
        compact_panels(0.01,0.01)
      end
      if doVertical
        hlinks = linkprop(h,{'YLim'});
        h(1).YLim = dep_lim;
        for ip = 2:npanels
          h(ip).YLabel.String = ''; 
          h(ip).YTickLabels = []; 
        end
      else
        hlinks = linkprop(h,{'XLim'});
        h(1).XLim = dep_lim;
      end
      
      set(gcf,'userdata',{'hlinks',hlinks})
      if nargout == 1
        varargout{1} = h;
      elseif nargout == 2
        varargout{1} = h;
        varargout{2} = hleg;
      elseif nargout == 2
        varargout{1} = h;
        varargout{2} = hleg;
        varargout{3} = hlinks;
      end
    end
    function varargout = plottimeseries(obj,varstrs_all,varargin)
       % Check additional input      
      have_options = 0;
      nargs = numel(varargin);      
      if nargs > 0, have_options = 1; args = varargin(:); end      
      while have_options
        l = 1;
        switch(lower(args{1}))
          case ''
          otherwise 
            warning(sprintf('Unknown argument %s.',args{1}))
        end
        args = args(l+1:end);  
        if isempty(args), break, end    
      end
           
      [nrows,ncols] = size(varstrs_all);    
      npanels = nrows*ncols;
      ip = 0;
      for irow = 1:nrows
        for icol = 1:ncols
          ip = ip + 1;
          h(irow,icol) = subplot(nrows,ncols,ip);
        end
      end
            
      hleg = gobjects(0);
      ip = 0;
      for irow = 1:nrows
        for icol = 1:ncols
          ip = ip + 1;
          isHoldOn = 0;
          %ip = sub2ind([ncols nrows],icol,irow);
          hca = h(irow,icol);
          % check if input demand som andditional input, e.g. n(1)
          varstrs = varstrs_all{ip};
          nvars = numel(varstrs);
          for ivar = 1:nvars
            if strfind(varstrs{ivar},'(')
              ind1 = strfind(varstrs{ivar},'(');
              ind2 = strfind(varstrs{ivar},')');            
              indstr = varstrs{ivar}(ind1+1:ind2-1);
              varstr =  varstrs{ivar}(1:ind1-1);
              var = obj.(varstr)(eval(indstr));
            else
              var = obj.(varstrs{ivar});
            end     
            if ivar == 2
              hold(hca,'on')
            end
            var = squeeze(mean(mean(var,1),2)); % mean over dimensions 1 and 2            
            displayname = [varstrs{ivar}];
            plot(hca,obj.twci,var,'DisplayName',displayname);              
            if not(isHoldOn)
              hold(hca,'on')
              isHoldOn = 1;
            end            
            hca.XLabel.String = ['t\omega_{ci}'];
          end
          hold(hca,'off')
          hca.XGrid = 'on';
          hca.YGrid = 'on';
          hleg(ip) = legend(hca,'location','eastoutside'); % this may make the panels of different widths, fix below
        end
      end
      drawnow;
      leftpos = 1;
      for ip = 1:npanels, panel_width(ip) = h(ip).Position(3); end      
      for ip = 1:npanels, h(ip).Position(3) = min(panel_width); end      
      titlestring = sprintf('x/d_i = [%.2f %.2f], z/d_i= [%.2f %.2f],',obj.xi(1),obj.xi(end),obj.zi(1),obj.zi(end));      
      h(1).Title.String = titlestring;
      compact_panels(0.01)
      hlinks = linkprop(h,{'XLim'});
      h(1).XLim = obj.twci([1 end]);
      set(gcf,'userdata',{'hlinks',hlinks})
      if nargout == 1
        varargout{1} = h;
      elseif nargout == 2
        varargout{1} = h;
        varargout{2} = hleg;
      elseif nargout == 2
        varargout{1} = h;
        varargout{2} = hleg;
        varargout{3} = hlinks;
      end
    end
              
    % Data analysis routines, time derivatives, interpolation, etc.
    % Interpolate fields    
    function varargout = interpline(obj,x0,z0,field,ds,varargin)
      % FIELDLINE Integrates field in direction of flow vector      
      %   interpline(obj,xlim,zlim,field,ds,varargin)
      %   Example: 
      %   x0 = 150;
      %   z0 = 9.6;
      %   ds = 0.01;      
      %   [linearclength,linex,liney,linez,linefx,linefy,linefz] =
      %   obj.interpline(xlim,zlim,field,ds,varargin);
       
      
      % Load field and make unit vectors
%       Fx = obj.xlim(x).zlim(z).([field 'x']);
%       Fy = obj.xlim(x).zlim(z).([field 'y']);
%       Fz = obj.xlim(x).zlim(z).([field 'z']);
%       fabs = sqrt(Fx.^2 + Fy.^2 + Fz.^2);
%       fx = Fx./fabs;
%       fy = Fy./fabs;
%       fz = Fz./fabs;
            
      % Plot progress (for debugging)
      doPlot = 0; 

      % Set default values if not given in input
      if isempty(ds); dx = 0.01; end
      
      % Define starting point
      arcline = 0;
      xline = x0;  
      yline = 0;
      zline = z0;
      
      % Interploate the magnetic field to starting point       
      fxline = obj.interp(xline(end),zline(end),obj.twpe,[field 'x']);
      fyline = obj.interp(xline(end),zline(end),obj.twpe,[field 'y']);
      fzline = obj.interp(xline(end),zline(end),obj.twpe,[field 'z']);

      if doPlot % Plot starting point
        plot(xline(end),zline(end),'.')
        drawnow;
        hold(gca,'on')
      end

      istep = 0;
      while 1
        istep = istep + 1;
        % Calculate step, magnetic field direction times the predefined stepsize
        xstep = fxline(istep)*ds;
        ystep = fyline(istep)*ds;
        zstep = fzline(istep)*ds;    
        % Advance line
        xline(istep+1) = xline(istep) + xstep;
        yline(istep+1) = yline(istep) + ystep;
        zline(istep+1) = zline(istep) + zstep;
        arcline(istep+1) = arcline(istep) + sqrt(xstep^2+ystep^2+zstep^2);

        % Break if we end up outside of box
%         if xline(end)<min([pic.xi(1) x(1)]) || ...
%            xline(end)>max([pic.xi(end) xlim(end)]) ||...
%            zline(end)<min([pic.zi(1) z(1)]) || ...
%            zline(end)>max([pic.zi(end) zlim(end)])
          if xline(end)<min(pic.xi(1)) || ...
           xline(end)>max(pic.xi(end)) ||...
           zline(end)<min(pic.zi(1)) || ...
           zline(end)>max(pic.zi(end))
          xline(end) = [];
          yline(end) = [];
          zline(end) = [];
          arcline(end) = [];
          break;
        end
        % Get field at new point    
        fxline(istep+1) = obj.interp(xline(end),zline(end),obj.twpe,[field 'x']);
        fyline(istep+1) = obj.interp(xline(end),zline(end),obj.twpe,[field 'y']);
        fzline(istep+1) = obj.interp(xline(end),zline(end),obj.twpe,[field 'z']);

        if doPlot
          plot(xline(end),zline(end),'.')
          drawnow;
        end

        if doArclengthLimit 
          if arcline(end)>nsteps
            break;
          end
        elseif istep > (nsteps-1)
          break;
        end

      end
      if doPlot
        hold(gca,'off')
      end


      varargout{1} = arcline;
      varargout{2} = xline;
      varargout{3} = yline;
      varargout{4} = zline;
      varargout{5} = fxline;
      varargout{6} = fyline;
      varargout{7} = fzline;     
    end
    function out = interpfield(obj,x,z,t,field,varargin)
      % Interpolates fields to given x,z,t
      %   out = interpfield(obj,x,z,t,field,varargin)
      %
      % x
       
      method = 'my_linear';
      method = 'linear';
      method = 'spline';
      
      nt = numel(t);
      sizex = size(x,1);
      sizez = size(z,1);
      idim_tx = sizex==nt;
      idim_tz = sizez==nt;
      nx = sizex(find(not(idim_tx)));
      nz = sizex(find(not(idim_tz)));
      
      
      
      nPoints = numel(x);      
      var = nan(size(x));      
      
      method = 2;
      switch method
        case 1
          for it = 1:nt
            for ix = 1:nx
              xmin = min(x(ix,it));
              xmax = max(x(ix,it));
              zmin = min(z(ix,it));
              zmax = max(z(ix,it));

              if find(t(it)==obj.twci)
                tmppic = obj.xlim([xmin xmax]).zlim(z(ix,it),'closest',2).twcilim(t(it),'exact');
                %data = obj.interpfield3(tmpx,tmpz,tmpt,tmpvar,x(it,:),z(it,:),t(it));
                tmpvar = tmppic.(field);
                tmpx   = tmppic.xi;
                tmpz   = tmppic.zi;
                new_data = interp1(tmpz,tmpvar',z(ix,it));
                var(ix,:,it) = new_data;
              else
                tmppic = obj.xlim([xmin xmax]).zlim([zmin zmax]).twcilim(t(it),'closest',2);


              end
            end
          end
        case 2
          for itime = 1:nt
            
            xmin = min(x(:,itime))-1;
            xmax = max(x(:,itime))+1;
            zmin = min(z(:,itime))-1;
            zmax = max(z(:,itime))+1;

            if find(t(itime)==obj.twci)
              tmppic = obj.xlim([xmin xmax]).zlim([zmin zmax]).twcilim(t(itime),'exact');
              %data = obj.interpfield3(tmpx,tmpz,tmpt,tmpvar,x(it,:),z(it,:),t(it));
              tmpvar = tmppic.(field);
              tmpx   = tmppic.xi;
              tmpz   = tmppic.zi;
              Vq = interp2(tmpx,tmpz,tmpvar',x(:,itime),z(:,itime));
              var(:,itime) = Vq;
            end
          end
      end          
      out = var;
    end
    function out = interp(obj,x,z,t,field,varargin)
      % Interpolates fields to given x,z,t
      %   out = interp(obj,x,z,t,field,varargin)
      %
      
      % Check input
      
      % Check if inteprolation is needed
      
      % Load bounding data
      
      % Interpolate field to a any number of point (x,z,t)
      %
      % PIC.INTERP_EB3 - Interpolate for a number of given points.
      %
      % To be implemented:
      %  - interpolation for several - SEEMS TO BE DONE
      %  - shape preserving interpolation
      %  - different interpolation types for temporal and spatial
      %    dimensions, particularly important for temporal dimension where
      %    the time steps are quite large
      
      method = 'linear';
      nBounding = 1;
      if strcmp(method,'linear')
        nClosest = 2;
      end
      
      nt = numel(t);
      if nt == 1
        nP = numel(x);
        t = repmat(t,nP,1);
      else nP = nt;
      end
%       sizex = size(x);
%       sizez = size(z);
%       idim_tx = sizex==nt;
%       idim_tz = sizez==nt;
%       nx = sizex(find(not(idim_tx)));
%       nz = sizex(find(not(idim_tz)));
      
      
      
      nPoints = nP;      
      var = nan(nPoints,1);
      
      for iP = 1:nPoints
        %iP
        % Make sure we get the pic object we need.
        %tmppic = get_pic_obj(obj);
        try
        if t(iP) == obj.twci(end) || t(iP) == obj.twci(1) % do 2D interpolation
          tmppic = obj.xlim(x(iP),'closest',nClosest).zlim(z(iP),'closest',nClosest).twcilim(t(iP));           
          tmpt =  tmppic.twci;
          tmpx =  tmppic.xi;
          tmpz =  tmppic.zi;          
          tmpvar = tmppic.get_exp(field); 
          [X,Z] = meshgrid(tmpx,tmpz);
          var(iP) = interp2(X,Z,tmpvar,x(iP),z(iP),method);
        elseif any(abs(obj.xi - x(iP)) < 1e-5) % end up very close to grid point, do 3 closest
          tmppic = obj.xlim(x(iP),'closest',3).zlim(z(iP),'closest',nClosest).twcilim(t(iP),'bounding',nBounding);
          tmpt =  tmppic.twci;
          tmpx =  tmppic.xi;
          tmpz =  tmppic.zi;          
          tmpvar = tmppic.get_exp(field); 
          [X,Z,T] = meshgrid(tmpx,tmpz,tmpt);
          var(iP) = interp3(X,Z,T,permute(tmpvar,[2 1 3]),x(iP),z(iP),t(iP),method);
         elseif any(abs(obj.zi - z(iP)) < 1e-5) % end up very close to grid point, do 3 closest
          tmppic = obj.xlim(x(iP),'closest',nClosest).zlim(z(iP),'closest',3).twcilim(t(iP),'bounding',nBounding);
          tmpt =  tmppic.twci;
          tmpx =  tmppic.xi;
          tmpz =  tmppic.zi;          
          tmpvar = tmppic.get_exp(field); 
          [X,Z,T] = meshgrid(tmpx,tmpz,tmpt);
          var(iP) = interp3(X,Z,T,permute(tmpvar,[2 1 3]),x(iP),z(iP),t(iP),method);
        else
          tmppic = obj.xlim(x(iP),'closest',nClosest).zlim(z(iP),'closest',nClosest).twcilim(t(iP),'bounding',nBounding);  
          tmpt =  tmppic.twci;
          tmpx =  tmppic.xi;
          tmpz =  tmppic.zi;
          tmpvar = tmppic.get_exp(field);
          [X,Z,T] = meshgrid(tmpx,tmpz,tmpt);
          var(iP) = interp3(X,Z,T,permute(tmpvar,[2 1 3]),x(iP),z(iP),t(iP),method);
        end
        % Check so no value ends up o the edges of, then the interpolation
        % doesn't work.
        
        % Interpolate to particle position
        % Vq = interp3(V,Xq,Yq,Zq) assumes X=1:N, Y=1:M, Z=1:P where [M,N,P]=SIZE(V).      

        %[X,Z,T] = meshgrid(tmpx,tmpz,tmpt);  
%         try
%         var(iP) = interp3(X,Z,T,permute(tmpvar,[2 1 3]),x(iP),z(iP),t(iP),method);
        catch
          1;
        end
        if isnan(var(iP))
          1;
        end
      end
      out = var;
      function out = get_pic_obj(pic)
        
      end
    end
    function out = get_points(obj,x,z,t,range,field,varargin)
      % Pic outs values in given spatial range around a given point
      %   out = get_points(obj,x,z,t,field,range,varargin)
      %
      nt = numel(t);
      if nt == 1
        nP = numel(x);
        t = repmat(t,nP,1);
      else
        nP = nt;
      end
      nPoints = nP;      
      var = nan(nPoints,1);
      
      for iP = 1:nPoints
        tmppic = obj.xlim(x(iP)+range).zlim(z(iP)+range).twcilim(t(iP));
        tmpt =  tmppic.twci;
        tmpx =  tmppic.xi;
        tmpz =  tmppic.zi;          
        tmpvar = tmppic.get_exp(field);         
        var(iP) = mean(mean(tmpvar));   
      end
      out = var;
    end
    function [Ex,Ey,Ez,Bx,By,Bz] = interp_EB(obj,x,z,t,varargin)
      % Interpolate field to a given point (x,z,t)
      % Used in paarticular for integrating particle trajectories. Proceed
      % with care if doing changes.
      %
      % To be implemented:
      %  - interpolation for several timesteps
      %  - shape preserving interpolation
      %  - different interpolation types for temporal and spatial
      %    dimensions, particularly important for temporal dimension where
      %    the time steps are quite large
            
      method = 'spline';
      %method = 'linear';
      
      if numel(varargin) > 0
        have_options = 1;
        args = varargin;
      end
      while have_options
        l = 1;
        switch lower(args{1})
          case 'spline'
            method = 'spline';
          case 'linear'
            method = 'linear';            
        end    
        args = args(l+1:end);
        if isempty(args), break, end 
      end
      
      switch method
        case 'spline'
          nBoundingT = 1; % spline for two points gives linear interpolation
          nClosestXZ = 5;
        case 'linear'
          nBoundingT = 1;
          nClosestXZ = 4;
      end
      
      
      nPoints = numel(t); 
      
      for iP = 1:nPoints
        if any(obj.twci==t(iP)) % exact match in time, do only 2D interp
          tmppic = obj.xlim(x(iP),'closest',nClosestXZ).zlim(z(iP),'closest',nClosestXZ).twcilim(t(iP),'exact');

          tmpt =  tmppic.twci;
          tmpx =  tmppic.xi;
          tmpz =  tmppic.zi;
          tmpEx = tmppic.Ex; %tmpEx(:,:,1) = smooth2(tmpEx(:,:,1),2,2); tmpEx(:,:,2) = smooth2(tmpEx(:,:,2),2,2);
          tmpEy = tmppic.Ey; 
          tmpEz = tmppic.Ez;
          tmpBx = tmppic.Bx;
          tmpBy = tmppic.By;
          tmpBz = tmppic.Bz;

          % Interpolate to particle position
          % Vq = interp3(V,Xq,Yq,Zq) assumes X=1:N, Y=1:M, Z=1:P where [M,N,P]=SIZE(V).      
          % The variables has their individual grids, so interpolate from these
          % individual grids to the common grid (Ey, and moments).
          [X_EX,Z_EX] = meshgrid(tmppic.xivar.Ex,tmppic.zivar.Ex);
          [X_EY,Z_EY] = meshgrid(tmppic.xivar.Ey,tmppic.zivar.Ey);
          [X_EZ,Z_EZ] = meshgrid(tmppic.xivar.Ez,tmppic.zivar.Ez);
          [X_BX,Z_BX] = meshgrid(tmppic.xivar.Bx,tmppic.zivar.Bx);
          [X_BY,Z_BY] = meshgrid(tmppic.xivar.By,tmppic.zivar.By);
          [X_BZ,Z_BZ] = meshgrid(tmppic.xivar.Bz,tmppic.zivar.Bz);

          try
          Ex(iP) = interp2(X_EX,Z_EX,permute(tmpEx,[2 1]),x(iP),z(iP),method);
          Ey(iP) = interp2(X_EY,Z_EY,permute(tmpEy,[2 1]),x(iP),z(iP),method);
          Ez(iP) = interp2(X_EZ,Z_EZ,permute(tmpEz,[2 1]),x(iP),z(iP),method);
          Bx(iP) = interp2(X_BX,Z_BX,permute(tmpBx,[2 1]),x(iP),z(iP),method);
          By(iP) = interp2(X_BY,Z_BY,permute(tmpBy,[2 1]),x(iP),z(iP),method);
          Bz(iP) = interp2(X_BZ,Z_BZ,permute(tmpBz,[2 1]),x(iP),z(iP),method);
          catch

            1;
          end          
        else
          tmppic = obj.xlim(x(iP),'closest',nClosestXZ).zlim(z(iP),'closest',nClosestXZ).twcilim(t(iP),'bounding',nBoundingT);          

          tmpt =  tmppic.twci;
          tmpx =  tmppic.xi;
          tmpz =  tmppic.zi;
          tmpEx = tmppic.Ex; %tmpEx(:,:,1) = smooth2(tmpEx(:,:,1),2,2); tmpEx(:,:,2) = smooth2(tmpEx(:,:,2),2,2);
          tmpEy = tmppic.Ey; 
          tmpEz = tmppic.Ez;
          tmpBx = tmppic.Bx;
          tmpBy = tmppic.By;
          tmpBz = tmppic.Bz;

          % Interpolate to particle position
          % Vq = interp3(V,Xq,Yq,Zq) assumes X=1:N, Y=1:M, Z=1:P where [M,N,P]=SIZE(V).      
          % The variables has their individual grids, so interpolate from these
          % individual grids to the common grid (Ey, and moments).
          [X_EX,Z_EX,T] = meshgrid(tmppic.xivar.Ex,tmppic.zivar.Ex,tmpt);
          [X_EY,Z_EY,T] = meshgrid(tmppic.xivar.Ey,tmppic.zivar.Ey,tmpt);
          [X_EZ,Z_EZ,T] = meshgrid(tmppic.xivar.Ez,tmppic.zivar.Ez,tmpt);
          [X_BX,Z_BX,T] = meshgrid(tmppic.xivar.Bx,tmppic.zivar.Bx,tmpt);
          [X_BY,Z_BY,T] = meshgrid(tmppic.xivar.By,tmppic.zivar.By,tmpt);
          [X_BZ,Z_BZ,T] = meshgrid(tmppic.xivar.Bz,tmppic.zivar.Bz,tmpt);

          try
          Ex(iP) = interp3(X_EX,Z_EX,T,permute(tmpEx,[2 1 3]),x(iP),z(iP),t(iP),method);
          Ey(iP) = interp3(X_EY,Z_EY,T,permute(tmpEy,[2 1 3]),x(iP),z(iP),t(iP),method);
          Ez(iP) = interp3(X_EZ,Z_EZ,T,permute(tmpEz,[2 1 3]),x(iP),z(iP),t(iP),method);
          Bx(iP) = interp3(X_BX,Z_BX,T,permute(tmpBx,[2 1 3]),x(iP),z(iP),t(iP),method);
          By(iP) = interp3(X_BY,Z_BY,T,permute(tmpBy,[2 1 3]),x(iP),z(iP),t(iP),method);
          Bz(iP) = interp3(X_BZ,Z_BZ,T,permute(tmpBz,[2 1 3]),x(iP),z(iP),t(iP),method);
          catch

            1;
          end          
        %plot(tmppic.xivar.Ez,tmpEz(:,:,2),'o',x(iP),Ez(iP),'x')
      %disp(sprintf('intEx = %g, meanEx = %g',Ex,mean(tmpEx(:))))
      % For debugging purposes
%       doPlot = 1;
%       if doPlot
%         figure(31)
%         if mod(nClosest,2) == 0 % even number          
%           hca = subplot(1,2,2);
%         else
%           hca = subplot(1,2,1);
%         end
%         scale = 2000;
%         scatter3(hca,X(:),Z(:),T(:),abs(tmpEx(:))*scale,abs(tmpEx(:)))
%         for ip = 1:numel(X)
%           text(hca,X(ip),Z(ip),T(ip),sprintf('%.3f',tmpEx(ip)))
%         end
%         hb = colorbar('peer',hca);
%         %plot3(hca,X(:),Z(:),T(:),'k.',x,z,t,'ro')
%         hold(hca,'on')
%         plot3(hca,x,z,t,'k+')
%         scatter3(hca,x,z,t,abs(Ex)*scale,abs(Ex))
%         hold(hca,'off')
%         hca.XLabel.String = 'x';
%         hca.YLabel.String = 'z';
%         hca.ZLabel.String = 't';
%         %pause
%       end
        end
      end
    end
    function [Ex,Ey,Ez,Bx,By,Bz] = interp_EB3(obj,x,z,t)
      % Interpolate field to a any number of points (x,z,t)
      %
      % PIC.INTERP_EB3 - Interpolate for a number of given points.
      %
      % To be implemented:
      %  - interpolation for several timesteps
      %  - shape preserving interpolation
      %  - different interpolation types for temporal and spatial
      %    dimensions, particularly important for temporal dimension where
      %    the time steps are quite large
      
      method = 'linear';
      if strcmp(method,'linear')
        % There is a problem when using 2, because some grid point are
        % offset by 0.5 from the "basic grid". It's slowed down considerably
        % though. 
        nClosestXZ = 3;
        % For time, 'closest' can become a problem with nonuniform spacing.
        % But, bounding also always returns 3 indices, which is
        % unnecessary. Changed bounding .
        nBoundingT = 1; 
        
      end
      
      nPoints = numel(t);
      
      Ex = nan(nPoints,1);
      Ey = nan(nPoints,1);
      Ez = nan(nPoints,1);
      Bx = nan(nPoints,1);
      By = nan(nPoints,1);
      Bz = nan(nPoints,1);
      
      for iP = 1:nPoints
      
        tmppic = obj.xlim(x(iP),'closest',nClosestXZ).zlim(z(iP),'closest',nClosestXZ).twcilim(t(iP),'bounding',nBoundingT);  

        tmpt =  tmppic.twci;
        tmpx =  tmppic.xi;
        tmpz =  tmppic.zi;
        tmpEx = tmppic.Ex;
        tmpEy = tmppic.Ey;
        tmpEz = tmppic.Ez;
        tmpBx = tmppic.Bx;
        tmpBy = tmppic.By;
        tmpBz = tmppic.Bz;
      
        % Interpolate to particle position
        % Vq = interp3(V,Xq,Yq,Zq) assumes X=1:N, Y=1:M, Z=1:P where [M,N,P]=SIZE(V).      

        %[X,Z,T] = meshgrid(tmpx,tmpz,tmpt);
        [X_EX,Z_EX,T] = meshgrid(tmppic.xivar.Ex,tmppic.zivar.Ex,tmpt);
        [X_EY,Z_EY,T] = meshgrid(tmppic.xivar.Ey,tmppic.zivar.Ey,tmpt);
        [X_EZ,Z_EZ,T] = meshgrid(tmppic.xivar.Ez,tmppic.zivar.Ez,tmpt);
        [X_BX,Z_BX,T] = meshgrid(tmppic.xivar.Bx,tmppic.zivar.Bx,tmpt);
        [X_BY,Z_BY,T] = meshgrid(tmppic.xivar.By,tmppic.zivar.By,tmpt);
        [X_BZ,Z_BZ,T] = meshgrid(tmppic.xivar.Bz,tmppic.zivar.Bz,tmpt);
              
        Ex(iP) = interp3(X_EX,Z_EX,T,permute(tmpEx,[2 1 3]),x(iP),z(iP),t(iP),method);
        Ey(iP) = interp3(X_EY,Z_EY,T,permute(tmpEy,[2 1 3]),x(iP),z(iP),t(iP),method);
        Ez(iP) = interp3(X_EZ,Z_EZ,T,permute(tmpEz,[2 1 3]),x(iP),z(iP),t(iP),method);
        Bx(iP) = interp3(X_BX,Z_BX,T,permute(tmpBx,[2 1 3]),x(iP),z(iP),t(iP),method);
        By(iP) = interp3(X_BY,Z_BY,T,permute(tmpBy,[2 1 3]),x(iP),z(iP),t(iP),method);
        Bz(iP) = interp3(X_BZ,Z_BZ,T,permute(tmpBz,[2 1 3]),x(iP),z(iP),t(iP),method); 
        if isnan(Ez(iP))
          disp(sprintf('ip = ',iP))
          disp(sprintf('xp = %g, xgrid = %g, %g, %g',x(iP),tmpx(1),tmpx(2),tmpx(3)))          
          disp(sprintf('zp = %g, zgrid = %g, %g, %g',z(iP),tmpz(1),tmpz(2),tmpz(3)))
          disp(sprintf('tp = %g, tgrid = %g, %g, ',t(iP),tmpt(1),tmpt(2)))
          1;
        end
      end
    end
    function [Ex,Ey,Ez,Bx,By,Bz] = interp_EB2(obj,x,z,t)
      % Interpolate field to a given point (x,z,t)
      %
      % PIC.INTERP_EB2 - Different interpolation for time and space 
      %
      % To be implemented:
      %  - interpolation for several timesteps
      %  - shape preserving interpolation
      %  - different interpolation types for temporal and spatial
      %    dimensions, particularly important for temporal dimension where
      %    the time steps are quite large
   
      nClosestXZ = 2;
      nClosestT = 5;            
      
      tmppic = obj.xlim(x,'closest',nClosestXZ).zlim(z,'closest',nClosestXZ).twcilim(t,'closest',nClosestT);
      
      tmpt =  tmppic.twci;
      tmpx =  tmppic.xi;
      tmpz =  tmppic.zi;
      tmpEx = tmppic.Ex;
      tmpEy = tmppic.Ey;
      tmpEz = tmppic.Ez;
      tmpBx = tmppic.Bx;
      tmpBy = tmppic.By;
      tmpBz = tmppic.Bz;
      
      % Interpolate to particle position
      % Vq = interp3(V,Xq,Yq,Zq) assumes X=1:N, Y=1:M, Z=1:P where [M,N,P]=SIZE(V).      

%      for iT = 1:nClosestT
        
      [X,Z,T] = meshgrid(tmpx,tmpz,tmpt);
      Ex = interp3(X,Z,T,permute(tmpEx,[2 1 3]),x,z,t,method);
      Ey = interp3(X,Z,T,permute(tmpEy,[2 1 3]),x,z,t,method);
      Ez = interp3(X,Z,T,permute(tmpEz,[2 1 3]),x,z,t,method);
      Bx = interp3(X,Z,T,permute(tmpBx,[2 1 3]),x,z,t,method);
      By = interp3(X,Z,T,permute(tmpBy,[2 1 3]),x,z,t,method);
      Bz = interp3(X,Z,T,permute(tmpBz,[2 1 3]),x,z,t,method);
      
   
    end
    function [vx,vy,vz] = interp_v(obj,x,z,t,iSpecies)
      % Interpolate field to a given point (x,z,t)
      %   [vx,vy,vz] = interp_v(obj,x,z,t,iSpecies)
      %
      % To be implemented:
      %  - interpolation for several timesteps
      %  - shape preserving interpolation
      %  - different interpolation types for temporal and spatial
      %    dimensions, particularly important for temporal dimension where
      %    the time steps are quite large
      
      method = 'linear';
      if strcmp(method,'linear')
        nClosest = 2;
      end
      
      nPoints = numel(t); 
      
      tmppic = obj.xlim(x,'closest',nClosest).zlim(z,'closest',nClosest).twcilim(t,'closest',nClosest);  
      
      tmpt =  tmppic.twci;
      tmpx =  tmppic.xi;
      tmpz =  tmppic.zi;
      tmpVx = tmppic.vx(iSpecies);
      tmpVy = tmppic.vy(iSpecies);
      tmpVz = tmppic.vz(iSpecies);

      % Interpolate to particle position
      % Vq = interp3(V,Xq,Yq,Zq) assumes X=1:N, Y=1:M, Z=1:P where [M,N,P]=SIZE(V).      

      [X,Z,T] = meshgrid(tmpx,tmpz,tmpt);
      vx = interp3(X,Z,T,permute(tmpVx,[2 1 3]),x,z,t,method);
      yy = interp3(X,Z,T,permute(tmpVy,[2 1 3]),x,z,t,method);
      zz = interp3(X,Z,T,permute(tmpVz,[2 1 3]),x,z,t,method);
      
    end
    function out = integrate_trajectory_old(obj,r0,v0,tspan,m,q,varargin)
      % out = integrate_trajectory(r0,v0,tspan,m,q,varargin)
      % tspan = [tstart tstop] - back or forward, if tstart > tstop, integrating is done backward in time
      %     or  [tstop_back tstart tstop_forw] -  integration is done forward and backward      
      % Additional input arguments varargin are passed on to as options in 
      % odeset. See 'help odeset', but commonly used are:
      %   'RelTol'
      %   'AbsTol'
      
      doPrintInfo = 0;
      
      if numel(tspan) == 2 % [tstart tstop]
        tstart = tspan(1);
        tstop_all = tspan(2);        
      elseif numel(tspan) == 3 % [tstop_back tstart tstop_forw]
        tstart = tspan(2);
        tstop_all = tspan([1 3]);
      end
      
      x_sol = [];
      for tstop = tstop_all
        % Print information
        if doPrintInfo
          if tstart < tstop, disp(['Integrating trajectory forward in time.'])
          else, disp(['Integrating trajectory backward in time.'])
          end
        end
        ttot = tic;
        x_init = [r0, v0]; % di, vA
        disp(sprintf('tstart = %5.2f, tstop = %5.2f, [x0,y0,z0] = [%5.1f, %5.1f, %5.1f], [vx0,vy0,vz0] = [%5.2f, %5.2f, %5.2f]',...
          tstart,tstop,x_init(1),x_init(2),x_init(3),x_init(4),x_init(5),x_init(6)))

        % Integrate trajectory
        %options = odeset();
        %options = odeset('AbsTol',1e-14,'Events',@exitBox);
        options = odeset('Events',@exitBox,varargin{:});
        %options = odeset('AbsTol',1e-7,'AbsTol',1e-9,'Events',@exitBox);
        %options = odeset('RelTol',1e-6);
        EoM = @(ttt,xxx) eom_pic(ttt,xxx,obj,m,q); 

        [t,x_sol_tmp] = ode45(EoM,[tstart tstop],x_init,options);%,options); % 
        x_sol_tmp(:,7) = t; % x_sol = (x,y,z,vx,vy,vz,t)

        x_sol = [x_sol; x_sol_tmp];
        
        doPlot = 0; % diagnostics
        if doPlot
          %%
          h = setup_subplots(2,2);

          hca = h(1);
          plot3(hca,x_sol(:,1),x_sol(:,2),x_sol(:,3),...
                    x_sol(1,1),x_sol(1,2),x_sol(1,3),'g*',...
                    x_sol(end,1),x_sol(end,2),x_sol(end,3),'r*')
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'y';        
          hca.ZLabel.String = 'z';
          hca.XGrid = 'on';
          hca.YGrid = 'on';
          hca.ZGrid = 'on';

          hca = h(2);
          plot(hca,x_sol(:,4),x_sol(:,5),...
                    x_sol(1,4),x_sol(1,5),'g*',...
                    x_sol(end,4),x_sol(end,5),'r*')
          hca.XLabel.String = 'vx';
          hca.YLabel.String = 'vy';
          hca.XGrid = 'on';
          hca.YGrid = 'on';

          hca = h(3);
          plot(hca,x_sol(:,4),x_sol(:,6),...
                    x_sol(1,4),x_sol(1,6),'g*',...
                    x_sol(end,4),x_sol(end,6),'r*')
          hca.XLabel.String = 'vx';
          hca.YLabel.String = 'vz';
          hca.XGrid = 'on';
          hca.YGrid = 'on';

          hca = h(4);
          plot(hca,x_sol(:,5),x_sol(:,6),...
                    x_sol(1,5),x_sol(1,6),'g*',...
                    x_sol(end,5),x_sol(end,6),'r*')
          hca.XLabel.String = 'vy';
          hca.YLabel.String = 'vz';
          hca.XGrid = 'on';
          hca.YGrid = 'on';
          drawnow
        end
        tt = toc(ttot);
        %out = x_sol;
      end
      % sort by time
      [~,I] = sort(x_sol(:,7)); % to implement: remove duplicates here
      out.t = x_sol(I,7);
      out.x = x_sol(I,1);
      out.y = x_sol(I,2);
      out.z = x_sol(I,3);
      out.vx = x_sol(I,4);
      out.vy = x_sol(I,5);
      out.vz = x_sol(I,6);   
      out.t0 = tstart;
      out.x0 = r0(1);
      out.y0 = r0(2);
      out.z0 = r0(3);
      out.vx0 = v0(1);
      out.vy0 = v0(2);
      out.vz0 = v0(3);
      out.options = options;      
    end    
    function out = integrate_trajectory(obj,r0,v0,tspan,m,q,varargin)
      % out = integrate_trajectory(r0,v0,tspan,m,q,varargin)
      % tspan = [tstart tstop] - back or forward, if tstart > tstop, integrating is done backward in time
      %     or  [tstop_back tstart tstop_forw] -  integration is done forward and backward      
      % Additional input arguments varargin are passed on to as options in 
      % odeset. See 'help odeset', but commonly used are:
      %   'RelTol'
      %   'AbsTol'
      
      % Default values
      doPrintInfo = 0;            
      method = 'spline';
      %method = 'linear';
      
%       if numel(varargin) > 0
%         have_options = 1;
%         args = varargin;
%       end
%       while have_options
%         l = 1;
%         switch lower(args{1})
%           case 'spline'
%             method = 'spline';
%           case 'linear'
%             method = 'linear';            
%         end    
%         args = args(l+1:end);
%         if isempty(args), break, end 
%       end
      disp(sprintf('Interpolation method: %s', method))
      if numel(tspan) == 2 % [tstart tstop]
        tstart = tspan(1);
        tstop_all = tspan(2);        
      elseif numel(tspan) == 3 % [tstop_back tstart tstop_forw]
        tstart = tspan(2);
        tstop_all = tspan([1 3]);
      end
      
      x_sol = [];
      x_sol_all = [];
      
      for tstop = tstop_all
        % Print information
        if doPrintInfo
          if tstart < tstop, disp(['Integrating trajectory forward in time.'])
          else, disp(['Integrating trajectory backward in time.'])
          end
        end
        ttot = tic;
        x_init = [r0, v0]; % di, vA
        disp(sprintf('tstart = %5.2f, tstop = %5.2f, [x0,y0,z0] = [%5.1f, %5.1f, %5.1f], [vx0,vy0,vz0] = [%5.2f, %5.2f, %5.2f]',...
          tstart,tstop,x_init(1),x_init(2),x_init(3),x_init(4),x_init(5),x_init(6)))

        % Integrate trajectory
        %options = odeset();
        %options = odeset('AbsTol',1e-14,'Events',@exitBox);
        %options = odeset('Events',@exitBox,varargin{:},'Events', @myEventCrossEquator);
        %options = odeset('Events',@exitBox,varargin{:},'Events', @(tt,xx) myEventCrossEquator(tt,xx,0));
        options = odeset('Events',@exitBox,varargin{:},'Events', @(tt,xx) myEventBoxEdge(tt,xx,obj.xi(([1 end]))));
        %options = odeset('Events',@exitBox,varargin{:},'Events', @(tt,xx) myEventBoxEdge(tt,xx,[190 210]));
        %options = odeset('AbsTol',1e-7,'AbsTol',1e-9,'Events',@exitBox);
        %options = odeset('RelTol',1e-6);
        EoM = @(ttt,xxx) eom_pic_internal(ttt,xxx,obj,m,q); 

        [t,x_sol_tmp] = ode45(EoM,[tstart tstop],x_init,options);%,options); % 
        x_sol_tmp(:,7) = t; % x_sol = (x,y,z,vx,vy,vz,t)

        x_sol = [x_sol; x_sol_tmp]; % 1st points are duplicates
        
        doPlot = 0; % diagnostics
        if doPlot
          %%
          h = setup_subplots(2,2);

          hca = h(1);
          plot3(hca,x_sol(:,1),x_sol(:,2),x_sol(:,3),...
                    x_sol(1,1),x_sol(1,2),x_sol(1,3),'g*',...
                    x_sol(end,1),x_sol(end,2),x_sol(end,3),'r*')
          hca.XLabel.String = 'x';
          hca.YLabel.String = 'y';        
          hca.ZLabel.String = 'z';
          hca.XGrid = 'on';
          hca.YGrid = 'on';
          hca.ZGrid = 'on';

          hca = h(2);
          plot(hca,x_sol(:,4),x_sol(:,5),...
                    x_sol(1,4),x_sol(1,5),'g*',...
                    x_sol(end,4),x_sol(end,5),'r*')
          hca.XLabel.String = 'vx';
          hca.YLabel.String = 'vy';
          hca.XGrid = 'on';
          hca.YGrid = 'on';

          hca = h(3);
          plot(hca,x_sol(:,4),x_sol(:,6),...
                    x_sol(1,4),x_sol(1,6),'g*',...
                    x_sol(end,4),x_sol(end,6),'r*')
          hca.XLabel.String = 'vx';
          hca.YLabel.String = 'vz';
          hca.XGrid = 'on';
          hca.YGrid = 'on';

          hca = h(4);
          plot(hca,x_sol(:,5),x_sol(:,6),...
                    x_sol(1,5),x_sol(1,6),'g*',...
                    x_sol(end,5),x_sol(end,6),'r*')
          hca.XLabel.String = 'vy';
          hca.YLabel.String = 'vz';
          hca.XGrid = 'on';
          hca.YGrid = 'on';
          drawnow
        end
        tt = toc(ttot);
        %out = x_sol;
      end
      
      % Sort by time
      [~,isort] = sort(x_sol_all(:,7));
      x_sol_all = x_sol_all(isort,:);
      
      % Find duplicate indices from interpolated values in solver
      idup = find(diff(x_sol_all(:,7))==0);
      x_sol_all(idup,:) = [];
      
      % sort by time
      [~,I] = sort(x_sol(:,7));                  
      out.t = x_sol(I,7);
      x_sol_all_ = interp1(x_sol_all(:,7),x_sol_all(:,1:6),x_sol(I,7));
      out.x = x_sol(I,1);
      out.y = x_sol(I,2);
      out.z = x_sol(I,3);
      out.vx = x_sol(I,4);
      out.vy = x_sol(I,5);
      out.vz = x_sol(I,6);   
      out.Ex = x_sol_all_(:,1);
      out.Ey = x_sol_all_(:,2);
      out.Ez = x_sol_all_(:,3);
      out.Bx = x_sol_all_(:,4);
      out.By = x_sol_all_(:,5);
      out.Bz = x_sol_all_(:,6);
      out.t0 = tstart;
      out.x0 = r0(1);
      out.y0 = r0(2);
      out.z0 = r0(3);
      out.vx0 = v0(1);
      out.vy0 = v0(2);
      out.vz0 = v0(3);
      out.options = options;   
      
      function  x_res = eom_pic_internal(t,x_vect,pic,m,q)
        x = x_vect(1);
        y = x_vect(2);
        z = x_vect(3);
        vx = x_vect(4);
        vy = x_vect(5);
        vz = x_vect(6);

        if isnan(x)
          1;
        end
        %disp(sprintf('t = %g, x = %g, y = %g, z = %g',t,x,y,z))

        %method = 'spline';        

        [Ex,Ey,Ez,Bx,By,Bz] = pic.interp_EB(x,z,t,method);
        
        %disp(sprintf('%.3f %.3f %.3f, %.3f %.3f %.3f, %.3f %.3f %.3f, %.3f, %.3f %.3f',Ex,Ey,Ez,Bx,By,Bz,x_vect(1),x_vect(2),x_vect(3),x_vect(4),x_vect(5),x_vect(6)))        
        it = size(x_sol_all,1);
        x_sol_all(it+1,1) = Ex;
        x_sol_all(it+1,2) = Ey;
        x_sol_all(it+1,3) = Ez;
        x_sol_all(it+1,4) = Bx;
        x_sol_all(it+1,5) = By;
        x_sol_all(it+1,6) = Bz;
        x_sol_all(it+1,7) = t;

        % Equations to be solved
        x_res = zeros(6,1);
        x_res(1) = vx; % dx/dt = vx;
        x_res(2) = vy; % dy/dt = vy;
        x_res(3) = vz; % dz/dt = vz;
        x_res(4) = (q/m)*(Ex + vy*Bz - vz*By);
        x_res(5) = (q/m)*(Ey + vz*Bx - vx*Bz);
        x_res(6) = (q/m)*(Ez + vx*By - vy*Bx);                                              

      end      
      function [value, isterminal, direction] = myEventBoxEdge(t, x_vect,boxedge)
        % integration is terminated when value changes sign
        % for this setup, value is initially negative
        value      = (x_vect(1)-boxedge(1))*(x_vect(1)-boxedge(2));        
        isterminal = 1;   % Stop the integration
        direction  = 0;
      end
      function [value, isterminal, direction] = myEventCrossEquator(t, x_vect,z)
        value      = x_vect(3) - z;
        isterminal = 1;   % Stop the integration
        direction  = 0;
      end
    end    
    function out = integrate_trajectory_constant_EB(obj,r0,v0,tstart,tstop,m,q)
      % out = integrate_trajectory(r0,v0,tstart,tstop,m,q)
      % if tstart > tstop, integrating is done backward in time
      
      T = tstop-tstart;
      
      if T > 0
        doForward = 1;
        %integration_string = 'forward';
        disp(['Integrating trajectory forward in time.'])
      else
        doForward = 0;
        %integration_string = 'backward';
        disp(['Integrating trajectory backward in time.'])        
      end
            
      
      ttot = tic;
      x_init = [r0, v0]; % di, vA
      disp(sprintf('tstart = %5.2f, tstop = %5.2f, [x0,y0,z0] = [%5.1f, %5.1f, %5.1f], [vx0,vy0,vz0] = [%5.2f, %5.2f, %5.2f]',...
        tstart,tstop,x_init(1),x_init(2),x_init(3),x_init(4),x_init(5),x_init(6)))

      % Integrate trajectory
      options = odeset('AbsTol',1e-14);
      [XI,ZI] = meshgrid(obj.xi,obj.zi);
      fields.x = XI;
      fields.z = ZI;
      fields.Bx = obj.Bx;
      fields.By = obj.By;
      fields.Bz = obj.Bz;
      fields.Ex = smooth_data(obj.Ex,10);
      fields.Ey = smooth_data(obj.Ey,10);
      fields.Ez = smooth_data(obj.Ez,10);
      
      if doForward
        EoM = @(ttt,xxx) eom_pic(ttt,xxx,fields,m,q);
      else
        %EoM = @(ttt,xxx) eom_pic_back(ttt,xxx,obj,m,q);
        EoM = @(ttt,xxx) eom_pic(ttt,xxx,fields,m,q);
      end        
      [t,x_sol] = ode45(EoM,[tstart tstop],x_init,options);%,options); % 
      x_sol(:,7) = t; % x_sol = (x,y,z,vx,vy,vz,t)

      doPlot = 1;
      if doPlot
        %%
        h = setup_subplots(2,2);
        
        hca = h(1);
        plot3(hca,x_sol(:,1),x_sol(:,2),x_sol(:,3),...
                  x_sol(1,1),x_sol(1,2),x_sol(1,3),'g*',...
                  x_sol(end,1),x_sol(end,2),x_sol(end,3),'r*')
        hca.XLabel.String = 'x';
        hca.YLabel.String = 'y';        
        hca.ZLabel.String = 'z';
        hca.XGrid = 'on';
        hca.YGrid = 'on';
        hca.ZGrid = 'on';
        
        hca = h(2);
        plot(hca,x_sol(:,4),x_sol(:,5),...
                  x_sol(1,4),x_sol(1,5),'g*',...
                  x_sol(end,4),x_sol(end,5),'r*')
        hca.XLabel.String = 'vx';
        hca.YLabel.String = 'vy';
        hca.XGrid = 'on';
        hca.YGrid = 'on';
        
        hca = h(3);
        plot(hca,x_sol(:,4),x_sol(:,6),...
                  x_sol(1,4),x_sol(1,6),'g*',...
                  x_sol(end,4),x_sol(end,6),'r*')
        hca.XLabel.String = 'vx';
        hca.YLabel.String = 'vz';
        hca.XGrid = 'on';
        hca.YGrid = 'on';
        
        hca = h(4);
        plot(hca,x_sol(:,5),x_sol(:,6),...
                  x_sol(1,5),x_sol(1,6),'g*',...
                  x_sol(end,5),x_sol(end,6),'r*')
        hca.XLabel.String = 'vy';
        hca.YLabel.String = 'vz';
        hca.XGrid = 'on';
        hca.YGrid = 'on';
        drawnow
      end
      toc(ttot)
      %out = x_sol;
      out.t = x_sol(:,7);
      out.x = x_sol(:,1);
      out.y = x_sol(:,2);
      out.z = x_sol(:,3);
      out.vx = x_sol(:,4);
      out.vy = x_sol(:,5);
      out.vz = x_sol(:,6);
    end
    function out = separatrix_location(obj,varargin)
      % PIC.SEPARATRIX_LOCATION Calculates location of separatrix
      %   At the moment it's the northern separatrix.
      %
      % sep = separatrix_location(pic);
      %      
      % sep = 
      %
      % struct with fields:
      %
      %      twci: [1×58 double]
      %         x: [58×6400 double]
      %         z: [58×6400 double]
      %   xline_x: [1×58 double]
      %   xline_z: [1×58 double]
      %
      % pcolor(sep.x(1,:),sep.twci,sep.z)
      % shading flat
      % hold on
      % plot(sep.xline_x,sep.twci,'k')
            
      doNorth = 1;
      nTimes = obj.nt;      
      x = obj.xi;
      z = obj.zi;
      %nTimes
      % Calculate separatrix for each time step
      new_x = nan(obj.nx,obj.nt);
      new_z = nan(obj.nx,obj.nt);
      for it = 1:nTimes
        %it
        twci = obj.twci(it);
        %twci
        A_tmp = squeeze(obj.twcilim(twci).A); % get A for this time step
        if mean(A_tmp(:)) < 0 % fix because I had wrong sign of A before
          A_tmp = -A_tmp;
        end
        [saddle_locations,saddle_values] = saddle(A_tmp,'sort'); % get saddle locations
        AX = saddle_values(1); % A-value of outermost X line
        
        % Save location of X line
        xX(it) = x(saddle_locations(1,1));
        zX(it) = z(saddle_locations(1,2));
        
        % Since the X line is a saddle point, it's very sensitive to small 
        % changes in A. Therefore, refine the grid around the X line and
        % find the location again.
        doRefine = 1;
        if doRefine
          % Indices of coarse X line.
          ixX(it) = saddle_locations(1,1);
          izX(it) = saddle_locations(1,2);
          
          % Refind A on a 21 x 21 grid around this point
          ind_x_sub = ixX(it)+(-10:10);
          ind_z_sub = izX(it)+(-10:10);
          n_old = numel(ind_x_sub); % number of original cells
          n_new = 100; % new number of cells for refinement
          x_old = x(ind_x_sub); % x values of original grid
          z_old = z(ind_z_sub); % z values of original grid
          A_old = A_tmp(ind_x_sub,ind_z_sub); % A values of original grid
          [X_old,Z_old] = meshgrid(x_old,z_old);
          x_new = linspace(x_old(1),x_old(end),n_new);
          z_new = linspace(z_old(1),z_old(end),n_new);
          [X_new,Z_new] = meshgrid(x_new,z_new);          
          % Interpolate A from old grid to new refined grid
          A_new = interp2(X_old,Z_old,A_old,X_new,Z_new);
          %contour(A_new')
          %drawnow
          %pause
          % when refining A, the MinPeakProminence needs to be lowered
          n_ref = n_old/n_new;
          default_minpeakprominence = 1e-2;
          foundSaddle = 0;
          while not(foundSaddle)
            [saddle_locations,saddle_values] = saddle(A_new,'sort','minpeakprominence',n_ref*default_minpeakprominence);
            if isempty(saddle_values)
              default_minpeakprominence = default_minpeakprominence*1e-1;
              disp(sprintf('Changed MinPeakProminence from %g to %g',default_minpeakprominence*1e1,default_minpeakprominence))
            else
              foundSaddle = 1;
              break;
            end
          end
          xX(it) = x_new(saddle_locations(1,1)); % new X line locations
          zX(it) = z_new(saddle_locations(1,2));
        end        
        % Pick a number of A values around the X line to get the magnetic 
        % field line from. Because In aleays want to find the northern one.
        % But sometimes I if just pick one A value, it becomes southern, or
        % inside the outflow (and then goes from north to south). From this
        % set of lines, I then find the one I want.
        dAmult = [1 1 1 1 1 1 1]+[-1e-2 -1e-4 -1e-6 0 1e-6 1e-4 1e-2];
        S = contourcs(x,z,A_tmp',AX*dAmult); % get contour lines
        % find contour of top (largest mean(z))
        nS = numel(S); % number of contour lines
        dxS = nan(nS,1);
        zminS = nan(nS,1);
        zmeanS = nan(nS,1);
        xmeanS = nan(nS,1);
        % Loop through lines and save some values. For example, a northern
        % line has larger mean(z) (Y) than a southern one, or one that is 
        % inside the outflow.
        for iS = 1:nS 
          dxS(iS) = S(iS).X(end) - S(iS).X(1);
          zminS(iS) = min(S(iS).Y);
          zmeanS(iS) = mean(S(iS).Y);
          xmeanS(iS) = mean(S(iS).X);
        end
        keep_ind_0 = 1:nS;        
        keep_ind_1 = find(dxS > 0.9*max(dxS)); % longest extent in x direction
        keep_ind_2 = find(zmeanS > 0); % should be above z = 0
        keep_ind_3 = intersect(keep_ind_1,keep_ind_2); % find lines that satisfy both these first two conditions
        keep_ind_4 = find(zminS == min(zminS(keep_ind_3))); % now find the line that has the smallest z value (because maybe theres two northern lines)
        keep_ind = keep_ind_4;
        indS = keep_ind;
               
        % remove all values below the x-line (in z)
        try % This was for bugfixing I think
        xSep = S(indS).X;
        zSep = S(indS).Y;
        catch
          1;
        end
          
        %keep_ind = find(zSep>=zX(it));
        %xSep = xSep(keep_ind);
        %zSep = zSep(keep_ind);
        
        if 1 % Diagnostic plotting
          hca = subplot(2,1,1);
          plot(hca,S(indS).X,S(indS).Y,xSep,zSep)
          hca = subplot(2,1,2);
          plot(hca,S(indS).X,atand(zSep./zSep))
          
          pause(0.1)
        end
        
        % Interpolate the chosen line to the simulation x-grid
        new_x(:,it) = x;
        new_z(:,it) = interp1(xSep,zSep,x); % only need to do this one ? no      
      end
      out.twci = obj.twci;
      out.x = new_x;
      out.z = new_z;
      out.xline_x = xX;
      out.xline_z = zX;
    end
    function out = intEpar(obj)
      % int(Epar)dlpar
      
      % Outer edge of box is zero, this might result in discrete jumps
      % where they meet. Start at upper half edge of box.
      pic = obj.zlim([0 obj.zi(end)]);
      intEpar = zeros(pic.nx,pic.nz);
      Epar = pic.Epar;
      Bx = pic.Bx;
      By = pic.By;
      Bz = pic.Bz;
      Babs = sqrt(Bx.^2 + By.^2 + Bz.^2);
      bx = Bx./Babs;
      by = By./Babs;
      bz = Bz./Babs;
      dx = (pic.xi(2)-pic.xi(1));
      dz = (pic.zi(2)-pic.zi(1));
      ds = dz/2;
      
      for ix = 2:pic.nx  
        1;
        disp(sprintf('ix = %g/%g',ix,pic.nx))
        for iz = 1:pic.nz          
          try            
            out = interp(obj,x,z,t,field,varargin);
          % interpolate from current to previous x          
          bx_ = bx(ix,iz);%*sign(bx(ix,iz));
          by_ = by(ix,iz);%*sign(bx(ix,iz));
          bz_ = bz(ix,iz);%*sign(bx(ix,iz));          
          % make interpolation in z-plane
          znew = pic.zi(iz)+dz;
          if znew > pic.zi(end) || znew < pic.zi(1)
            continue
          end
          iz1 = find(pic.zi<znew,1,'last');
          iz2 = find(pic.zi>znew,1,'first');
          
          % distance to points
          dz1 = (pic.zi(iz1)-znew)/dz;
          dz2 = (pic.zi(iz2)-znew)/dz;
          % Values at the bounding pounts
%           Bx1 = Bx(ix,iz1); Bx2 = Bx(ix,iz2);
%           By1 = By(ix,iz1); By2 = By(ix,iz2);
%           Bz1 = Bz(ix,iz1); Bz2 = Bz(ix,iz2);
          Epar1 = Epar(ix,iz1); Epar2 = Epar(ix,iz2);
          intEpar1 = intEpar(ix,iz1); intEpar2 = intEpar(ix,iz2);
          
%           Bxprev = dz2 * Bx1 + dz1 * Bx2;
%           Bxprev = dz2 * By1 + dz1 * By2;
%           Bxprev = dz2 * Bz1 + dz1 * Bz2;
          Eparprev = dz2 * Epar1 + dz1 * Epar2;
          intEparprev = dz2 * intEpar1 + dz1 * intEpar2;
          intEpar(ix,iz) = intEparprev-Eparprev*ds;
          catch
            1;
          end
        end
      end
      out = intEpar;
      
    end
    
    % Get simulation meta data and parameters 
    function out = parse_namelist(obj)
      % PARSE_NAMELIST Find basic information of simulation.
      % textscan
      % regexp
      
      buffer = fileread(obj.namelist);
      % Simulation information
      out.mime = str2double(regexpi(buffer, '(?<=mime\s*=\s*)\d*', 'match'));
      out.wpewce = str2double(regexpi(buffer, '(?<=wpewce\s*=\s*)\d*', 'match'));
      out.c = str2double(regexpi(buffer, '(?<=c\s*=\s*)\d*', 'match'));     
      out.teti = str2double(regexpi(buffer, '(?<=theta\s*=\s*)(\d*\.\d*|\d*)', 'match'));
      % Plasma species information
      names = regexpi(buffer, '(?<=name\s*=\s*)''(\w*)''', 'match');
      for iname = 1:numel(names)
        names{iname} = strrep(names{iname},'''','');
      end
      out.name = names;
      mass_ = regexpi(buffer, '(?<=mass\s*=\s*)(\w*|\d*)', 'match');
      for imass = 1:numel(mass_)
      if strcmp(mass_{imass},'mime')
        out.mass(imass) = out.mime;
      else
        out.mass(imass) = str2double(mass_{imass});
      end
      end
        
      out.charge = str2double(regexpi(buffer, '(?<=charge\s*=\s*)(-\d*|\d*)', 'match'));
      out.particles_per_cell = str2double(regexpi(buffer, '(?<=particles_per_cell\s*=\s*)(-\d*|\d*)', 'match'));
      
      % expression for box size
      [~,gr] = grep('-s',regexpi(buffer,'\nbox_size\s*=\s*','match'),obj.namelist);
      str_box_size = gr.match{1};
      str_box_size = strrep(str_box_size,'m.','');
      str_box_size = strrep(str_box_size,'mime','out.mime');
      eval([str_box_size,';']);
      resx = str2double(regexpi(buffer, '(?<=resx\s*=\s*)\d*', 'match'));
      resy = str2double(regexpi(buffer, '(?<=resy\s*=\s*)\d*', 'match'));
      resz = str2double(regexpi(buffer, '(?<=resy\s*=\s*)\d*', 'match'));
      cell_length = [1./resx, 1./resy, 1./resz];
      grid_length = box_size;
      out.nx = box_size(1)/cell_length(1);
      out.ny = box_size(2)/cell_length(2);
      out.nz = box_size(3)/cell_length(3);
      out.xe = 0:cell_length(1):(box_size(1));
      out.ye = 0:cell_length(2):(box_size(2));
      out.ze = 0:cell_length(3):(box_size(3));
      
      % particle binning
      deposited_quantity = regexpi(buffer, '(?<=deposited_quantity\s*=\s*)("\w*")', 'match');      
      species = regexpi(buffer, '(?<=species\s*=\s*[)("\w*")', 'match');
       for ii = 1:numel(deposited_quantity)
        deposited_quantity{ii} = strrep(deposited_quantity{ii},'"','');
        species{ii} = strrep(species{ii},'"','');
        species{ii} = strrep(species{ii},'[','');
       end
      out.deposited_quantity = deposited_quantity;
      out.deposited_species = species;            
      
    end
    function out = get_software(obj)
      nAttr = numel(obj.info.Attributes);
      for iAttr = 1:nAttr
        %attributes{iAttr} = obj.info.Attributes(iAttr).Name;
        if strcmp(obj.info.Attributes(iAttr).Name,'software')
          out = obj.info.Attributes(iAttr).Value;
          return
        end
      end
    end
    function out = get_attributes(obj)
      %nAttr = numel(obj.info.Attributes);
      allAttr = {};
      nIter = obj.nt;
      for iIter = 1:nIter
        nAttr = numel(obj.info.Groups(1).Groups(iIter).Attributes);
        for iAttr = 1:nAttr
          allAttr{end+1} = obj.info.Groups(1).Groups(iIter).Attributes(iAttr).Name;
        end
          %attributes.(obj.info.Groups(1).Attributes(iAttr).Name) = obj.info.Groups(1).Attributes(iAttr).Value;        
        end
      out = unique(allAttr);
    end
    function out = get_timeline_attributes(obj,attr_str)
      fileInfo = obj.info;
      iGroup = find(contains({fileInfo.Groups.Name},'/data'));
      nIter = numel(fileInfo.Groups(iGroup).Groups); % number of iterations
      
      isub = 0;
      for iIter = obj.it
        isub = isub + 1;
        % /data/00000xxxxx/ 
        % redo to actually find the 
        %iIter
        iAtt = find(contains({fileInfo.Groups(iGroup).Groups(iIter).Attributes.Name},attr_str));
        if numel(iAtt)>1
          if iIter == 1
            warning(sprintf('Found more than one (partial) match to %s, choosing the first match %s.',attr_str,fileInfo.Groups(iGroup).Groups(iIter).Attributes(iAtt(1)).Name))
          end
          iAtt = iAtt(1);
        end
        if not(isempty(iAtt))
          datasize = fileInfo.Groups(iGroup).Groups(iIter).Attributes(iAtt).Dataspace.Size;          
          attr(isub,:) = fileInfo.Groups(iGroup).Groups(iIter).Attributes(iAtt).Value;
        else
          attr(isub,:) = NaN;
        end
        
      end
      out = attr;
    end
    function out = get_missing_attributes(obj,attr_str)
      fileInfo = obj.info;
      iGroup = find(contains({fileInfo.Groups.Name},'/data'));
      nIter = numel(fileInfo.Groups(iGroup).Groups); % number of iterations
      
      for iIter = 1:nIter
        % /data/00000xxxxx/ 
        % redo to actually find the 
        %iIter
        iAtt = find(contains({fileInfo.Groups(iGroup).Groups(iIter).Attributes.Name},attr_str));
        if isempty(iAtt)
          isMissing(iIter) = 1;
        else 
          isMissing(iIter) = 0;          
        end
        
      end
      out = find(isMissing);
    end
    function out = get_twpe(obj)
      fileInfo = obj.info;
      iGroup = find(contains({fileInfo.Groups.Name},'/data'));
      nIter = numel(fileInfo.Groups(iGroup).Groups); % number of iterations
      
      for iIter = 1:nIter
        % /data/00000xxxxx/ 
        % redo to actually find the 
        iAtt = find(contains({fileInfo.Groups(iGroup).Groups(iIter).Attributes.Name},'time'));
        time(iIter) = fileInfo.Groups(iGroup).Groups(iIter).Attributes(iAtt).Value;
      end
      precision = 10; % after I change the timestep for nobg, there were some very small number in the 9th digit
      time = fix(time*precision)/precision;
      out = time;
    end
    function out = get_iterations(obj)
      fileInfo = obj.info;
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
      %%
      % needs to be adapted for the species subgroups
      fileInfo = obj.info;
      % fields structure is the same for all times, but different for
      % micPIC and Smilei      
      if strcmp(obj.software,'micPIC') % plasma moments are in subgroup with datasets for each species
        fields = {fileInfo.Groups(1).Groups(1).Datasets.Name};
        split_str = cellfun(@(x) strsplit(x,'/'),{fileInfo.Groups(1).Groups(1).Groups.Name},'UniformOutput',false);
        for iout = 1:numel(split_str)
          fields{end+1} = split_str{iout}{end};
        end
      elseif strcmp(obj.software,'Smilei') % all basic data are in the same group, 2nd order moments are in different files
        all_fields = {fileInfo.Groups(1).Groups(1).Datasets.Name};
        namelist = parse_namelist(obj);        
        species = namelist.name;      
        ifield = 1;
        for iOut = 1:numel(all_fields)              
          if strfind(all_fields{iOut},'_') % is plasma field
            tokens = strsplit(all_fields{iOut},'_');
            % check if field already exist            
            tmp = cellfun(@(x)logical(strfind(x,tokens{1})),fields,'UniformOutput',false);            
            if isempty([tmp{:}])                       
              fields{ifield} = tokens{1};
            else
              continue;
            end
          else
            fields{ifield} = all_fields{iOut};
          end
          ifield = ifield + 1;
        end
        % load particle binning fields
        
      end
      out = fields;
    end
    function out = get_gridsize(obj)
      out = [numel(obj.xe) numel(obj.ye) numel(obj.ze)];
    end
    function out = get_mass(obj)
      fileInfo = obj.info;
      if strcmp(obj.software,'micPIC')
        iGroup = find(contains({fileInfo.Groups.Name},'/data'));
        iAtt = find(contains({fileInfo.Groups(iGroup).Groups(1).Groups(1).Datasets(1).Attributes.Name},'mass'));
        nSpecies = numel(fileInfo.Groups(iGroup).Groups(1).Groups(1).Datasets);
        for iSpecies = 1:nSpecies
          mass(iSpecies) = fileInfo.Groups(iGroup).Groups(1).Groups(1).Datasets(iSpecies).Attributes(iAtt).Value;
        end
        out = mass;
      else
        out = obj.mass;
      end
      out = out; % mass is normalized to mass of first species.
    end
    function out = get_charge(obj)
      fileInfo = obj.info;
      if strcmp(obj.software,'micPIC')
        iGroup = find(contains({fileInfo.Groups.Name},'/data'));
        iAtt = find(contains({fileInfo.Groups(iGroup).Groups(1).Groups(1).Datasets(1).Attributes.Name},'charge'));
        nSpecies = numel(fileInfo.Groups(iGroup).Groups(1).Groups(1).Datasets);
        for iSpecies = 1:nSpecies
        charge(iSpecies) = fileInfo.Groups(iGroup).Groups(1).Groups(1).Datasets(iSpecies).Attributes(iAtt).Value;
        end
        out = charge;
      else 
        out = obj.charge; % this was already loaded in using namelist
      end      
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
    function out = ny(obj)
      out = numel(obj.yi);
    end
    function out = nz(obj)
      out = numel(obj.zi);
    end
    function out = nt(obj)
      out = obj.length;
    end
    
    % Get fields
    % Convention is to build in the coordinate tranformation when loading
    % the Smilei data: x -> x, -z -> y, y -> z
    function out = A(obj)
      % pic.A Magnetic vector potential Ay.
      if any(contains(obj.fields,'A')) % stored as field
       	out = get_field(obj,'A');
      else % calculate it
         if strcmp(obj.software,'micPIC')
          Bx =  obj.Bx;
          Bz =  obj.Bz;
          A = calc_A(obj.xi,obj.zi,Bx,Bz);
        elseif strcmp(obj.software,'Smilei')
          Bx =  obj.Bx; 
          Bz =  obj.Bz;% coordinate tranformation built-in into .Bz and .By
          A = calc_A(obj.xi,obj.zi,Bx,Bz);
        end
        out = A;
      end
      % nested function
      function out = calc_A(x,z,bx,bz)
        % Grid
        dx = x(2)-x(1);
        dz = z(2)-z(1);
        nz = numel(z);
        nx = numel(x);


        ntimes = size(bx,3);
        A = bx*0;
        for itime = 1:ntimes
          bx_tmp = squeeze(bx(:,:,itime));
          bz_tmp = squeeze(bz(:,:,itime));
          A_tmp = zeros(nx,nz);
          % Dont put zero right at the edge 
          ixm = 10;    
          % Advance up
          A_tmp(ixm,:) = cumsum(bx_tmp(ixm,:),2)*dz;
          % Advance to the right
          A_tmp((ixm+1):end,:) = repmat(A_tmp(ixm,:),nx-ixm,1) + cumsum(-bz_tmp((ixm+1):end,:),1)*dx;
          % Advance to the left
          A_tmp((ixm-1):-1:1,:) = repmat(A_tmp(ixm,:),ixm-1,1) + cumsum(-bz_tmp((ixm-1):-1:1,:),1)*dx;
          A(:,:,itime) = A_tmp;
        end
        out = A;
      end
    end
    function out = Ax(obj)
      % Assume Ax is defined by By. By = - dAz/dx + dAx/dz, with Az = 0
      
      By = obj.By;
      %By = smooth2(By,100);
      
      A = calc_Ax(obj.xi,obj.zi,By);
        
      out = A;      
      % nested function
      function out = calc_Ax(x,z,by)
        % Grid
        dx = x(2)-x(1);
        dz = z(2)-z(1);
        nz = numel(z);
        nx = numel(x);


        ntimes = size(by,3);
        A = by*0;
        for itime = 1:ntimes
          by_tmp = squeeze(by(:,:,itime));          
          % Advance up
          A_tmp = cumsum(by_tmp,2)*dz;                    
          A(:,:,itime) = A_tmp;
        end
        out = A;
      end
    end
    function out = Az(obj)
      % Assume Az is defined by By. By = dAz/dx - dAx/dz, with Ax = 0
      
      By =  obj.By;
      
      A = calc_Az(obj.xi,obj.zi,By);
        
      out = A;      
      % nested function
      function out = calc_Az(x,z,by)
        % Grid
        dx = x(2)-x(1);
        dz = z(2)-z(1);
        nz = numel(z);
        nx = numel(x);


        ntimes = size(by,3);
        A = by*0;
        for itime = 1:ntimes
          by_tmp = squeeze(by(:,:,itime));          
          % Advance up
          A_tmp = -cumsum(by_tmp,1)*dx;
          A(:,:,itime) = A_tmp;
        end
        out = A;
      end
    end
    function out = Bx(obj)
      if strcmp(obj.software,'micPIC')
        out = get_field(obj,'bx')*obj.wpewce;
      elseif strcmp(obj.software,'Smilei')
        out = get_field(obj,'Bx')*obj.wpewce;
      end
    end
    function out = By(obj)
      if strcmp(obj.software,'micPIC')
        out = get_field(obj,'by')*obj.wpewce;
      elseif strcmp(obj.software,'Smilei')
        out = get_field(obj,'By')*obj.wpewce;
      end
    end
    function out = Bz(obj)
      if strcmp(obj.software,'micPIC')
        out = get_field(obj,'bz')*obj.wpewce;
      elseif strcmp(obj.software,'Smilei')
        out = get_field(obj,'Bz')*obj.wpewce;
      end
    end
    function out = Babs(obj)
      Bx = obj.Bx;
      By = obj.By;
      Bz = obj.Bz;
      out = sqrt(Bx.^2 + By.^2 + Bz.^2);     
    end
    function out = Ex(obj)
      if strcmp(obj.software,'micPIC')
        out = get_field(obj,'ex')*sqrt(obj.mime)*obj.wpewce^2;
      elseif strcmp(obj.software,'Smilei')
        out = get_field(obj,'Ex')*sqrt(obj.mime)*obj.wpewce^2;
      end      
    end
    function out = Ey(obj)
      if strcmp(obj.software,'micPIC')
        out = get_field(obj,'ey')*sqrt(obj.mime)*obj.wpewce^2;
      elseif strcmp(obj.software,'Smilei')
        out = get_field(obj,'Ey')*sqrt(obj.mime)*obj.wpewce^2;
      end      
    end
    function out = Ez(obj)
      if strcmp(obj.software,'micPIC')
        out = get_field(obj,'ez')*sqrt(obj.mime)*obj.wpewce^2;
      elseif strcmp(obj.software,'Smilei')
        out = get_field(obj,'Ez')*sqrt(obj.mime)*obj.wpewce^2;
      end     
    end
    function out = PB(obj)
      % Magnetic field pressure
      %   out = PB(obj,value)      
      Bx = obj.Bx;
      By = obj.By;
      Bz = obj.Bz;
      out = 0.5*(Bx.^2 + By.^2 + Bz.^2);
    end
    function out = pB(obj)
      % Magnetic field pressure
      %   out = pB(obj,value)            
      out = obj.PB;
    end
    % Density
    function out = n(obj,species)
      % pic.N Density of selected populations defined by indices.
      %   iSpecies = [2 4]; n = pic.N(iSpecies);
      %   n = pic.N(1);
      %   n = pic.N([1 3]);
      % See also PIC.NI, PIC.NE
      
      % Check that only a single species is given (unique charge)
      % Different masses are ok (for example protons and oxygen)
      if not(numel(unique(obj.charge(species))) == 1)
        error('Selected species have different mass and/or charge. This is not supported.')
      end
      nSpecies = numel(species);
      
      if strcmp(obj.software,'micPIC')
        dfac = obj.get_dfac;
        n = zeros(obj.nx,obj.nz,obj.nt);
        for iSpecies = species
          n = n + obj.get_field(sprintf('dns/%.0f',iSpecies))*dfac(iSpecies);
        end        
      elseif strcmp(obj.software,'Smilei')
        n = zeros(obj.nx,obj.nz,obj.nt);
        for iSpecies = species
          pop_str = obj.species{iSpecies};
          n = n + obj.charge(iSpecies)*obj.get_field(['Rho_' pop_str]); % Smilei stores charge density, so to get density we multiply with -1
        end    
      end      
      out = n;
    end
    function out = ne(obj)
      % pic.NE Total electron density.      
      %   ne = pic.NE;      
      % See also PIC.N, PIC.NI
      iSpecies = find(obj.get_charge == -1); % negatively charge particles are electrons
      out = obj.n(iSpecies);
    end
    function out = ni(obj)
      % pic.NI Total ion density.
      %   ni = pic.NI;      
      % See also PIC.N, PIC.NE
      iSpecies = find(obj.get_charge == 1); % negatively charge particles are electrons
      out = obj.n(iSpecies);
    end
    % Flux
    function out = jx(obj,species)
      % pic.JX Particle flux jx = n*vx of selected populations.      
      %     iSpecies = [2 4]; jx = pic.JX(iSpecies);
      %     jx = pic.JX(1);
      %     jx = pic.JX([1 3]);
      %   To get partial current, you need to multiply with charge (1,-1).      
      %   Jx = pic.JIX - pic.JEX;
      % See also: PIC.JY, PIC.JZ, PIC.JEX, PIC.JEY, PIC.JEZ, PIC.JIX,
      %   PIC.JIY, PIC.JIZ
      
      % Check that only a single species of a given charge is given
      if not(numel(unique(obj.charge(species))) == 1)
        error('Selected species have different charge. This is not supported.')
      end
      nSpecies = numel(species);
      if strcmp(obj.software,'micPIC')
        dfac = obj.get_dfac;
        var = zeros(obj.nx,obj.nz,obj.nt);
        for iSpecies = species
          var = var + obj.get_field(sprintf('vxs/%.0f',iSpecies))*dfac(iSpecies)*obj.wpewce*sqrt(obj.mime);
        end
      elseif strcmp(obj.software,'Smilei')
        var = zeros(obj.nx,obj.nz,obj.nt);
        for iSpecies = species
          pop_str = obj.species{iSpecies};
          var = var + obj.charge(iSpecies)*obj.get_field(['Jx_' pop_str])*obj.wpewce*sqrt(obj.mime); % normalization ???
        end
      end
      out = var;
    end
    function out = jy(obj,species)
      % Get jy
      
      % Check that only a single species of a given charge is given
      if not(numel(unique(obj.charge(species))) == 1)
        error('Selected species have different charge. This is not supported.')
      end
      nSpecies = numel(species);
      if strcmp(obj.software,'micPIC')
        dfac = obj.get_dfac;
        var = zeros(obj.nx,obj.nz,obj.nt);
        for iSpecies = species
          var = var + obj.get_field(sprintf('vys/%.0f',iSpecies))*dfac(iSpecies)*obj.wpewce*sqrt(obj.mime);
        end
      elseif strcmp(obj.software,'Smilei')
        var = zeros(obj.nx,obj.nz,obj.nt);
        for iSpecies = species
          pop_str = obj.species{iSpecies};
          var = var + obj.charge(iSpecies)*obj.get_field(['Jy_' pop_str])*obj.wpewce*sqrt(obj.mime); % normalization ???
        end
      end
      out = var;    
    end
    function out = jz(obj,species)
      % Get jz
      
      % Check that only a single species of a given charge is given
      if not(numel(unique(obj.charge(species))) == 1)
        error('Selected species have different charge. This is not supported.')
      end
      nSpecies = numel(species);
      if strcmp(obj.software,'micPIC')
        dfac = obj.get_dfac;
        var = zeros(obj.nx,obj.nz,obj.nt);
        for iSpecies = species
          var = var + obj.get_field(sprintf('vzs/%.0f',iSpecies))*dfac(iSpecies)*obj.wpewce*sqrt(obj.mime);
        end
      elseif strcmp(obj.software,'Smilei')
        var = zeros(obj.nx,obj.nz,obj.nt);
        for iSpecies = species
          pop_str = obj.species{iSpecies};
          var = var + obj.charge(iSpecies)*obj.get_field(['Jz_' pop_str])*obj.wpewce*sqrt(obj.mime); % normalization ???
        end
      end
      out = var;  
    end
    function out = jex(obj)
      % Get electron flux, x
      iSpecies = find(obj.get_charge == -1); % negatively charge particles are electrons
      out = obj.jx(iSpecies);
    end
    function out = jey(obj)
      % Get electron flux, y
      iSpecies = find(obj.get_charge == -1); % negatively charge particles are electrons
      out = obj.jy(iSpecies);
    end
    function out = jez(obj)
      % Get electron flux, z
      iSpecies = find(obj.get_charge == -1); % negatively charge particles are electrons
      out = obj.jz(iSpecies);
    end
    function out = jix(obj)
      iSpecies = find(obj.get_charge == 1); % negatively charge particles are electrons
      out = obj.jx(iSpecies);
    end
    function out = jiy(obj)
      iSpecies = find(obj.get_charge == 1); % negatively charge particles are electrons
      out = obj.jy(iSpecies);
    end
    function out = jiz(obj)
      iSpecies = find(obj.get_charge == 1); % negatively charge particles are electrons
      out = obj.jz(iSpecies);
    end
    % Velocity
    function out = vx(obj,species)
      % pic.VX Velocity vx of selected populations.
      %     iSpecies = [2 4]; vx = pic.VX(iSpecies);
      %     vx = pic.VX(1);
      %     vx = pic.VX([1 3]);     
      n = obj.n(species);
      jx = obj.jx(species);      
      out = jx./n;      
    end
    function out = vy(obj,species)
      % Get velocity vy          
      n = obj.n(species);
      jy = obj.jy(species);      
      out = jy./n;      
    end
    function out = vz(obj,species)
      % Get velocity vz            
      n = obj.n(species);
      jz = obj.jz(species);      
      out = jz./n;      
    end
    function out = vabs(obj,species)
      vx = obj.vx(species);
      vy = obj.vy(species);
      vz = obj.vz(species);
      out = sqrt(vx.^2 + vy.^2 + vz.^2);
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
    function out = vix(obj)
      % Get ion velocity, x
      out = obj.jix./obj.ni;
    end
    function out = viy(obj)
      % Get ion velocity, y
      out = obj.jiy./obj.ni;
    end
    function out = viz(obj)
      % Get ion velocity, z
      out = obj.jiz./obj.ni;
    end   
    % Current
    function out = Jx(obj)
      out = obj.jix - obj.jex;
    end
    function out = Jy(obj)
      out = obj.jiy - obj.jey;
    end
    function out = Jz(obj)
      out = obj.jiz - obj.jez;
    end
    % Stress tensor
    function out = vexx(obj)
      iSpecies = find(obj.get_charge == -1); % negatively charge particles are electrons
      dfac = obj.get_dfac;      
      var = zeros([obj.get_gridsize,1]);
      for iComp = 1:numel(iSpecies)
        dataset = sprintf('vxx/%.0f',iSpecies(iComp));
        var_tmp = get_field(obj,dataset);
        var = var + var_tmp*dfac(iComp)*mass(iSpecies)*wpewce^2;
      end
      out = var;
      out = [];
    end
    function out = vxx(obj,species)
      % Stress tensor components sum(mvv)
      %   out = vxx(obj,species)
      if strcmp(obj.software,'micPIC')
        dfac = obj.get_dfac;         
        dataset = sprintf('vxx/%.0f',species);
        out = obj.mass(species)*obj.wpewce^2*get_field(obj,dataset)*dfac(species); 
      elseif strcmp(obj.software,'Smilei')
        out = obj.mass(species)*obj.wpewce^2*obj.get_binned_quantity('weight_vx_px',species);        
      end
    end
    function out = vxy(obj,species)
      % Stress tensor components sum(mvv)
      %   out = vxy(obj,species)
      if strcmp(obj.software,'micPIC')
        dfac = obj.get_dfac;         
        dataset = sprintf('vxy/%.0f',species);
        out = obj.mass(species)*obj.wpewce^2*get_field(obj,dataset)*dfac(species); 
      elseif strcmp(obj.software,'Smilei')
        out = obj.get_binned_quantity('weight_vx_py',species);        
      end
    end
    function out = vxz(obj,species)
      % Stress tensor components sum(mvv)
      %   out = vxz(obj,species)
      if strcmp(obj.software,'micPIC')
        dfac = obj.get_dfac;         
        dataset = sprintf('vxz/%.0f',species);
        out = obj.mass(species)*obj.wpewce^2*get_field(obj,dataset)*dfac(species); 
      elseif strcmp(obj.software,'Smilei')
        out = -obj.get_binned_quantity('weight_vx_pz',species); % added minus sign for coordinate transformation        
      end
    end
    function out = vyy(obj,species)
      % Stress tensor components sum(mvv)
      %   out = vyy(obj,species)
      
      if strcmp(obj.software,'micPIC')
        dfac = obj.get_dfac;
        dataset = sprintf('vyy/%.0f',species);
        out = obj.mass(species)*obj.wpewce^2*get_field(obj,dataset)*dfac(species);
      elseif strcmp(obj.software,'Smilei')
        out = obj.get_binned_quantity('weight_vy_py',species);
      end
    end
    function out = vyz(obj,species)
      % Stress tensor components sum(mvv)
      %   out = vyz(obj,species)
      if strcmp(obj.software,'micPIC')
        dfac = obj.get_dfac;
        dataset = sprintf('vyz/%.0f',species);
        out = obj.mass(species)*obj.wpewce^2*get_field(obj,dataset)*dfac(species);
      elseif strcmp(obj.software,'Smilei')
        out = -obj.get_binned_quantity('weight_vy_pz',species); % added minus sign for coordinate transformation        
      end
    end
    function out = vzz(obj,species)
      % Stress tensor components sum(mvv)
      %   out = vzz(obj,species)
      if strcmp(obj.software,'micPIC')
        dfac = obj.get_dfac;
        dataset = sprintf('vzz/%.0f',species);
        out = obj.mass(species)*obj.wpewce^2*get_field(obj,dataset)*dfac(species);
      elseif strcmp(obj.software,'Smilei')
        out = obj.mass(species)*obj.wpewce^2*obj.get_binned_quantity('weight_vz_pz',species);
      end      
    end
    function out = vv_diag(obj,species)
      %   out = vv_diag(obj,value)
%       dfac = obj.get_dfac;         
%       vxx = get_field(obj,sprintf('vxx/%.0f',species))*dfac(species);
%       vyy = get_field(obj,sprintf('vyy/%.0f',species))*dfac(species);
%       vzz = get_field(obj,sprintf('vzz/%.0f',species))*dfac(species);
%       out = obj.mass(species)*obj.wpewce^2*(vxx + vyy + vzz)/3;
      out = (obj.vxx(species) + obj.vyy(species) + obj.vzz(species))/3;
    end
    function out = pDxx(obj,species)
      % Get dynamical pressure mn(vx^2 + vy^2 + vz^2)/3
      % Kinetic energy is (3/2)*p_dyn (=mv^2/2)
      
      mass_sp = obj.mass; 
      if numel(unique(mass_sp(species))) > 1
        error('All species do not have the same mass.'); 
      else
        mass_sp = mass_sp(species(1));
      end
      charge = obj.get_charge; charge(species);
      if numel(unique(charge(species))) > 1
        error('All species do not have the same charge.');         
      end
      
      dfac = obj.get_dfac;
      % Initialize variables
      n = zeros(numel(obj.xi),numel(obj.zi),obj.length);
      vs = zeros(numel(obj.xi),numel(obj.zi),obj.length);
      
      % Sum over species
      for iSpecies = species                
        n = n + obj.get_field(sprintf('dns/%.0f',iSpecies))*dfac(iSpecies);  % density      
        vs = vs + obj.get_field(sprintf('vxs/%.0f',iSpecies))*dfac(iSpecies); % flux        
      end        
      % p_dyn =  mnvv
      out = mass_sp*obj.wpewce^2*(vs.*vs)./n/3;
    end
    function out = pDyy(obj,species)
      % Get dynamical pressure mn(vx^2 + vy^2 + vz^2)/3
      % Kinetic energy is (3/2)*p_dyn (=mv^2/2)
      
      mass_sp = obj.mass; 
      if numel(unique(mass_sp(species))) > 1
        error('All species do not have the same mass.'); 
      else
        mass_sp = mass_sp(species(1));
      end
      charge = obj.get_charge; charge(species);
      if numel(unique(charge(species))) > 1
        error('All species do not have the same charge.');         
      end
      
      dfac = obj.get_dfac;
      % Initialize variables
      n = zeros(numel(obj.xi),numel(obj.zi),obj.length);
      vs = zeros(numel(obj.xi),numel(obj.zi),obj.length);
      
      % Sum over species
      for iSpecies = species                
        n = n + obj.get_field(sprintf('dns/%.0f',iSpecies))*dfac(iSpecies);  % density      
        vs = vs + obj.get_field(sprintf('vxs/%.0f',iSpecies))*dfac(iSpecies); % flux        
      end        
      % p_dyn =  mnvv
      out = mass_sp*obj.wpewce^2*(vs.*vs)./n/3;
    end
    function out = pDzz(obj,species)
      % Get dynamical pressure mn(vx^2 + vy^2 + vz^2)/3
      % Kinetic energy is (3/2)*p_dyn (=mv^2/2)
      
      mass_sp = obj.mass; 
      if numel(unique(mass_sp(species))) > 1
        error('All species do not have the same mass.'); 
      else
        mass_sp = mass_sp(species(1));
      end
      charge = obj.get_charge; charge(species);
      if numel(unique(charge(species))) > 1
        error('All species do not have the same charge.');         
      end
      
      dfac = obj.get_dfac;
      % Initialize variables
      n = zeros(numel(obj.xi),numel(obj.zi),obj.length);
      vs = zeros(numel(obj.xi),numel(obj.zi),obj.length);
      
      % Sum over species
      for iSpecies = species                
        n = n + obj.get_field(sprintf('dns/%.0f',iSpecies))*dfac(iSpecies);  % density      
        vs = vs + obj.get_field(sprintf('vxs/%.0f',iSpecies))*dfac(iSpecies); % flux        
      end        
      % p_dyn =  mnvv
      out = mass_sp*obj.wpewce^2*(vs.*vs)./n/3;
    end
    function out = pD(obj,species)
      % Get dynamical pressure mn(vx^2 + vy^2 + vz^2)/3
      % Kinetic energy is (3/2)*p_dyn (=mv^2/2)
      
      mass_sp = obj.mass; 
      if numel(unique(mass_sp(species))) > 1
        error('All species do not have the same mass.'); 
      else
        mass_sp = mass_sp(species(1));
      end
      charge = obj.get_charge; charge(species);
      if numel(unique(charge(species))) > 1
        error('All species do not have the same charge.');         
      end
      
      dfac = obj.get_dfac;
      % Initialize variables
      n = zeros(numel(obj.xi),numel(obj.zi),obj.length);
      vxs = zeros(numel(obj.xi),numel(obj.zi),obj.length);
      vys = zeros(numel(obj.xi),numel(obj.zi),obj.length);
      vzs = zeros(numel(obj.xi),numel(obj.zi),obj.length);
      
      % Sum over species
      for iSpecies = species                
        n = n + obj.get_field(sprintf('dns/%.0f',iSpecies))*dfac(iSpecies);  % density      
        vxs = vxs + obj.get_field(sprintf('vxs/%.0f',iSpecies))*dfac(iSpecies); % flux
        vys = vys + obj.get_field(sprintf('vys/%.0f',iSpecies))*dfac(iSpecies); % flux
        vzs = vzs + obj.get_field(sprintf('vzs/%.0f',iSpecies))*dfac(iSpecies); % flux
      end        
      % p_dyn =  mnvv
      out = mass_sp*obj.wpewce^2*(vxs.*vxs + vys.*vys + vzs.*vzs)./n/3;
    end
    % Pressure, not implemented for Smilei
    function out = p12(obj,species,comp)
      % p12 = p12(obj,species,comp)
      % load general tensor components, 
      % comp is some of 'xx','xy','xz','yy','yz','zz'
      
      comp = sort(comp);
      
      nSpecies = numel(species);
      % check so that all species have the same mass and charge
      mass_sp = obj.mass; 
      if numel(unique(mass_sp(species))) > 1
        error('All species do not have the same mass.'); 
      else
        mass_sp = mass_sp(species(1));
      end
      charge = obj.get_charge; charge(species);
      if numel(unique(charge(species))) > 1
        error('All species do not have the same charge.');
      end
        
      if strcmp(obj.software,'micPIC')
        dfac = obj.get_dfac;

        % Initialize variables
        n = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        v1s = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        v2s = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        v12 = zeros(numel(obj.xi),numel(obj.zi),obj.length);

        % Sum over species
        for iSpecies = species
          n = n + obj.get_field(sprintf('dns/%.0f',iSpecies))*dfac(iSpecies);  % density
          v1s = v1s + obj.get_field(sprintf('v%ss/%.0f',comp(1),iSpecies))*dfac(iSpecies); % flux 
          if strcmp(comp(1),comp(2))
            v2s = v1s;
          else
            v2s = v2s + obj.get_field(sprintf('v%ss/%.0f',comp(2),iSpecies))*dfac(iSpecies); % flux
          end
          v12 = v12 + obj.get_field(sprintf('v%s/%.0f',comp,iSpecies))*dfac(iSpecies); % nvv
        end
        % p = P - mnvv
        out = mass_sp*obj.wpewce^2*(v12 - v1s.*v2s./n); % p12
      elseif strcmp(obj.software,'Smilei')
        % Initialize variables
        n = zeros(numel(obj.xi),numel(obj.zi),obj.length);      
        j1 = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        j2 = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        v12 = zeros(numel(obj.xi),numel(obj.zi),obj.length);      

        % Sum over species
        for iSpecies = species
          n = n + obj.n(iSpecies);  % density      
          eval(sprintf('j1 = j1 + obj.j%s(iSpecies);',comp(1))) % flux     
          if strcmp(comp(1),comp(2))
            j2 = j1;
          else
            eval(sprintf('j2 = j2 + obj.j%s(iSpecies);',comp(2))) % flux     
          end
          eval(sprintf('v12 = v12 + obj.v%s(iSpecies);',comp)) % nvv
        end        
        % p = P - mnvv
        out = mass_sp*obj.wpewce^2*(v12 - j1.*j2./n); % p12
        
      end
    end
    % To reduce lines of code, make these following call p12 when Smilei 
    % has been verified.
    function out = pxx(obj,species)
      % pxx = pxx(obj,species)
      
      nSpecies = numel(species);
      % check so that all species have the same mass and charge
      mass_sp = obj.mass; 
      if numel(unique(mass_sp(species))) > 1
        error('All species do not have the same mass.'); 
      else
        mass_sp = mass_sp(species(1));
      end
      charge = obj.get_charge; charge(species);
      if numel(unique(charge(species))) > 1
        error('All species do not have the same charge.');         
      end
        
      if strcmp(obj.software,'micPIC')
        dfac = obj.get_dfac;

        % Initialize variables
        n = zeros(numel(obj.xi),numel(obj.zi),obj.length);      
        vxs = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vxx = zeros(numel(obj.xi),numel(obj.zi),obj.length);      

        % Sum over species
        for iSpecies = species
          n = n + obj.get_field(sprintf('dns/%.0f',iSpecies))*dfac(iSpecies);  % density      
          vxs = vxs + obj.get_field(sprintf('vxs/%.0f',iSpecies))*dfac(iSpecies); % flux     
          vxx = vxx + obj.get_field(sprintf('vxx/%.0f',iSpecies))*dfac(iSpecies); % nvv
        end        
        % p = P - mnvv
        out = mass_sp*obj.wpewce^2*(vxx - vxs.*vxs./n); % pxx
      elseif strcmp(obj.software,'Smilei')
        % Initialize variables
        n = zeros(numel(obj.xi),numel(obj.zi),obj.length);      
        jx = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vxx = zeros(numel(obj.xi),numel(obj.zi),obj.length);      

        % Sum over species
        for iSpecies = species
          n = n + obj.n(iSpecies);  % density      
          jx = jx + obj.jx(iSpecies); % flux     
          vxx = vxx + obj.vxx(iSpecies); % nvv
        end        
        % p = P - mnvv
        out = mass_sp*obj.wpewce^2*(vxx - jx.*jx./n); % pxx
      end
    end
    function out = pxy(obj,species)
      % pxx = pxx(obj,species)
      
      nSpecies = numel(species);
      % check so that all species have the same mass and charge
      mass_sp = obj.mass; 
      if numel(unique(mass_sp(species))) > 1
        error('All species do not have the same mass.'); 
      else
        mass_sp = mass_sp(species(1));
      end
      charge = obj.get_charge; charge(species);
      if numel(unique(charge(species))) > 1
        error('All species do not have the same charge.');         
      end
        
      if strcmp(obj.software,'micPIC')
        dfac = obj.get_dfac;

        % Initialize variables
        n = zeros(numel(obj.xi),numel(obj.zi),obj.length);      
        vxs = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vys = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vxy = zeros(numel(obj.xi),numel(obj.zi),obj.length);      

        % Sum over species
        for iSpecies = species
          n = n + obj.get_field(sprintf('dns/%.0f',iSpecies))*dfac(iSpecies);  % density      
          vxs = vxs + obj.get_field(sprintf('vxs/%.0f',iSpecies))*dfac(iSpecies); % flux     
          vys = vys + obj.get_field(sprintf('vys/%.0f',iSpecies))*dfac(iSpecies); % flux     
          vxy = vxy + obj.get_field(sprintf('vxy/%.0f',iSpecies))*dfac(iSpecies); % nvv
        end        
        % p = P - mnvv
        out = mass_sp*obj.wpewce^2*(vxy - vxs.*vys./n); % pxy
      elseif strcmp(obj.software,'Smilei')
        % not implemented
      end
    end
    function out = pxz(obj,species)
      % pxx = pxx(obj,species)
      
      nSpecies = numel(species);
      % check so that all species have the same mass and charge
      mass_sp = obj.mass; 
      if numel(unique(mass_sp(species))) > 1
        error('All species do not have the same mass.'); 
      else
        mass_sp = mass_sp(species(1));
      end
      charge = obj.get_charge; charge(species);
      if numel(unique(charge(species))) > 1
        error('All species do not have the same charge.');         
      end
        
      if strcmp(obj.software,'micPIC')
        dfac = obj.get_dfac;

        % Initialize variables
        n = zeros(numel(obj.xi),numel(obj.zi),obj.length);      
        vxs = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vzs = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vxz = zeros(numel(obj.xi),numel(obj.zi),obj.length);      

        % Sum over species
        for iSpecies = species
          n = n + obj.get_field(sprintf('dns/%.0f',iSpecies))*dfac(iSpecies);  % density      
          vxs = vxs + obj.get_field(sprintf('vxs/%.0f',iSpecies))*dfac(iSpecies); % flux     
          vzs = vzs + obj.get_field(sprintf('vzs/%.0f',iSpecies))*dfac(iSpecies); % flux     
          vxz = vxz + obj.get_field(sprintf('vxz/%.0f',iSpecies))*dfac(iSpecies); % nvv
        end        
        % p = P - mnvv
        out = mass_sp*obj.wpewce^2*(vxz - vxs.*vzs./n); % pxy
      elseif strcmp(obj.software,'Smilei')
        % not implemented
      end
    end
    function out = pyz(obj,species)
      % pxx = pxx(obj,species)
      
      nSpecies = numel(species);
      % check so that all species have the same mass and charge
      mass_sp = obj.mass; 
      if numel(unique(mass_sp(species))) > 1
        error('All species do not have the same mass.'); 
      else
        mass_sp = mass_sp(species(1));
      end
      charge = obj.get_charge; charge(species);
      if numel(unique(charge(species))) > 1
        error('All species do not have the same charge.');         
      end
        
      if strcmp(obj.software,'micPIC')
        dfac = obj.get_dfac;

        % Initialize variables
        n = zeros(numel(obj.xi),numel(obj.zi),obj.length);      
        vys = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vzs = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vyz = zeros(numel(obj.xi),numel(obj.zi),obj.length);      

        % Sum over species
        for iSpecies = species
          n = n + obj.get_field(sprintf('dns/%.0f',iSpecies))*dfac(iSpecies);  % density      
          vys = vys + obj.get_field(sprintf('vys/%.0f',iSpecies))*dfac(iSpecies); % flux     
          vzs = vzs + obj.get_field(sprintf('vzs/%.0f',iSpecies))*dfac(iSpecies); % flux     
          vyz = vyz + obj.get_field(sprintf('vyz/%.0f',iSpecies))*dfac(iSpecies); % nvv
        end        
        % p = P - mnvv
        out = mass_sp*obj.wpewce^2*(vyz - vys.*vzs./n); % pyz
      elseif strcmp(obj.software,'Smilei')
        % not implemented
      end
    end
    function out = pyy(obj,species)
      % pxx = pxx(obj,species)
      
      nSpecies = numel(species);
      % check so that all species have the same mass and charge
      mass_sp = obj.mass; 
      if numel(unique(mass_sp(species))) > 1
        error('All species do not have the same mass.'); 
      else
        mass_sp = mass_sp(species(1));
      end
      charge = obj.get_charge; charge(species);
      if numel(unique(charge(species))) > 1
        error('All species do not have the same charge.');         
      end
        
      if strcmp(obj.software,'micPIC')
        dfac = obj.get_dfac;

        % Initialize variables
        n = zeros(numel(obj.xi),numel(obj.zi),obj.length);      
        vs = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vv = zeros(numel(obj.xi),numel(obj.zi),obj.length);      

        % Sum over species
        for iSpecies = species
          n = n + obj.get_field(sprintf('dns/%.0f',iSpecies))*dfac(iSpecies);  % density      
          vs = vs + obj.get_field(sprintf('vys/%.0f',iSpecies))*dfac(iSpecies); % flux     
          vv = vv + obj.get_field(sprintf('vyy/%.0f',iSpecies))*dfac(iSpecies); % nvv
        end        
        % p = P - mnvv
        out = mass_sp*obj.wpewce^2*(vv - vs.*vs./n); % pxx
      elseif strcmp(obj.software,'Smilei')
        % Initialize variables
        n = zeros(numel(obj.xi),numel(obj.zi),obj.length);      
        j = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vv = zeros(numel(obj.xi),numel(obj.zi),obj.length);      

        % Sum over species
        for iSpecies = species
          n = n + obj.n(iSpecies);  % density      
          j = j + obj.jy(iSpecies); % flux     
          vv = vv + obj.vyy(iSpecies); % nvv
        end        
        % p = P - mnvv
        out = mass_sp*obj.wpewce^2*(vv - j.*j./n); % pzz
      end
    end
    function out = pzz(obj,species)
      % pxx = pxx(obj,species)
      
      nSpecies = numel(species);
      % check so that all species have the same mass and charge
      mass_sp = obj.mass; 
      if numel(unique(mass_sp(species))) > 1
        error('All species do not have the same mass.'); 
      else
        mass_sp = mass_sp(species(1));
      end
      charge = obj.get_charge; charge(species);
      if numel(unique(charge(species))) > 1
        error('All species do not have the same charge.');         
      end
        
      if strcmp(obj.software,'micPIC')
        dfac = obj.get_dfac;

        % Initialize variables
        n = zeros(numel(obj.xi),numel(obj.zi),obj.length);      
        vxs = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vxx = zeros(numel(obj.xi),numel(obj.zi),obj.length);      

        % Sum over species
        for iSpecies = species
          n = n + obj.get_field(sprintf('dns/%.0f',iSpecies))*dfac(iSpecies);  % density      
          vxs = vxs + obj.get_field(sprintf('vzs/%.0f',iSpecies))*dfac(iSpecies); % flux     
          vxx = vxx + obj.get_field(sprintf('vzz/%.0f',iSpecies))*dfac(iSpecies); % nvv
        end        
        % p = P - mnvv
        out = mass_sp*obj.wpewce^2*(vxx - vxs.*vxs./n); % pxx
      elseif strcmp(obj.software,'Smilei')
        % Initialize variables
        n = zeros(numel(obj.xi),numel(obj.zi),obj.length);      
        j = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vv = zeros(numel(obj.xi),numel(obj.zi),obj.length);      

        % Sum over species
        for iSpecies = species
          n = n + obj.n(iSpecies);  % density      
          j = j + obj.jz(iSpecies); % flux     
          vv = vv + obj.vzz(iSpecies); % nvv
        end        
        % p = P - mnvv
        out = (vv - j.*j./n); % pzz
      end
    end
    function out = p(obj,species)
      % p = (pxx+pyy+pzz)/3
      out = (obj.pxx(species) + obj.pyy(species) + obj.pzz(species))/3;
    end
    function out = pi(obj,species)
      % p = (pxx+pyy+pzz)/3
      species = find(obj.get_charge == 1); % negatively charge particles are electrons      
      out = (obj.pxx(species) + obj.pyy(species) + obj.pzz(species))/3;
    end
    function out = pe(obj,species)
      % p = (pxx+pyy+pzz)/3
      species = find(obj.get_charge == -1); % negatively charge particles are electrons      
      out = (obj.pxx(species) + obj.pyy(species) + obj.pzz(species))/3;
    end
    function out = pixx(obj)
      % p = (pxx+pyy+pzz)/3
      species = find(obj.get_charge == 1); % negatively charge particles are electrons      
      out = obj.pxx(species);
    end
    function out = pixy(obj)
      % p = (pxx+pyy+pzz)/3
      species = find(obj.get_charge == 1); % negatively charge particles are electrons      
      out = obj.pxy(species);
    end
    function out = pixz(obj)
      % p = (pxx+pyy+pzz)/3
      species = find(obj.get_charge == 1); % negatively charge particles are electrons      
      out = obj.pxz(species);
    end
    function out = piyy(obj)
      % p = (pxx+pyy+pzz)/3
      species = find(obj.get_charge == 1); % negatively charge particles are electrons      
      out = obj.pyy(species);
    end
    function out = piyz(obj)
      % p = (pxx+pyy+pzz)/3
      species = find(obj.get_charge == 1); % negatively charge particles are electrons      
      out = obj.pyz(species);
    end
    function out = pizz(obj)
      % p = (pxx+pyy+pzz)/3
      species = find(obj.get_charge == 1); % negatively charge particles are electrons      
      out = obj.pzz(species);
    end
    function out = pexx(obj)
      % p = (pxx+pyy+pzz)/3
      species = find(obj.get_charge == -1); % negatively charge particles are electrons      
      out = obj.pxx(species);
    end
    function out = pexy(obj)
      % p = (pxx+pyy+pzz)/3
      species = find(obj.get_charge == -1); % negatively charge particles are electrons      
      out = obj.pxy(species);
    end
    function out = pexz(obj)
      % p = (pxx+pyy+pzz)/3
      species = find(obj.get_charge == -1); % negatively charge particles are electrons      
      out = obj.pxz(species);
    end
    function out = peyy(obj)
      % p = (pxx+pyy+pzz)/3
      species = find(obj.get_charge == -1); % negatively charge particles are electrons      
      out = obj.pyy(species);
    end
    function out = peyz(obj)
      % p = (pxx+pyy+pzz)/3
      species = find(obj.get_charge == -1); % negatively charge particles are electrons      
      out = obj.pyz(species);
    end
    function out = pezz(obj)
      % p = (pxx+pyy+pzz)/3
      species = find(obj.get_charge == -1); % negatively charge particles are electrons      
      out = obj.pzz(species);
    end
    % Temperature, not implemented for Smilei
    function out = t12(obj,species,comp)
      p = obj.p12(species,comp);
      n = obj.n(species);
      out = p./n;
    end
    function out = t(obj,species)
      n = obj.n(species);
      p = obj.p(species);
      out = p./n;
    end
    function out = te(obj)
      species = find(obj.get_charge == -1); % negatively charge particles are electrons      
      n = obj.n(species);
      p = obj.p(species);
      out = p./n;
    end
    function out = ti(obj)
      species = find(obj.get_charge == 1); % negatively charge particles are electrons      
      n = obj.n(species);
      p = obj.p(species);
      out = p./n;
    end
    function out = t_diag(obj,species)
      n = obj.n(species);
      pxx = obj.pxx(species);
      pyy = obj.pyy(species);
      pzz = obj.pzz(species);      
      t.xx = pxx./n;
      t.yy = pyy./n;
      t.zz = pzz./n;
      out = t;
    end
    function out = t_tens(obj,species)
      n = obj.n(species);
      pxx = obj.pxx(species);
      pxy = obj.pxy(species);
      pxz = obj.pxz(species);
      pyy = obj.pyy(species);
      pyz = obj.pyz(species);
      pzz = obj.pzz(species);
      t.xx = pxx./n;
      t.xy = pxy./n;
      t.xz = pxz./n;
      t.yy = pyy./n;
      t.yz = pyz./n;
      t.zz = pzz./n;
      out = t;
    end
    function out = t_fac(obj,species)
      % PIC.T_FAC Load t_tens, and rotate to field aligned coordinate system
      % tfac = PIC.T_FAC(species)
      % r1 = B/|B|;
      % r2 = r1 x [0 1 0] - in inflow, without guide field, this is then
      %                     close to z
      % r3 = r1 x r2
      
      % Temperature
      t = obj.t_tens(species);
      % Magnetic field
      Bx = obj.Bx;
      By = obj.By;
      Bz = obj.Bz;
      Babs = sqrt(Bx.^2 + By.^2 + Bz.^2);
      b.x =  Bx./Babs;
      b.y =  By./Babs;
      b.z =  Bz./Babs;
      % New coordinate system
      r1 = b; % magnetic field unit vector
      r2 = cross_product(r1.x,r1.y,r1.z,0,1,0);
      r2.abs = sqrt(r2.x.^2 + r2.y.^2 + r2.z.^2);
      r2.x = r2.x./r2.abs;
      r2.y = r2.y./r2.abs;
      r2.z = r2.z./r2.abs;
      r2.abs = sqrt(r2.x.^2 + r2.y.^2 + r2.z.^2);
      r2 = cross_product(r2.x,r2.y,r2.z,r1.x,r1.y,r1.z);
      r3 = cross_product(r1.x,r1.y,r1.z,r2.x,r2.y,r2.z);
      r3.abs = sqrt(r3.x.^2 + r3.y.^2 + r3.z.^2);
      
      % Rotate tensor
      % To get fac, we dont really need the entire tensor, do we? The 
      % diagonal should be enough?
      t_fac = rotate_tens(t,r1,r2,r3);
      t_perp = 0.5*(t_fac.yy + t_fac.zz);
      t_par = t_fac.xx;  
      t_scal = (t_fac.xx + t_fac.yy + t_fac.zz)/3;
      out.fac = t_fac;
      out.perp = t_perp;
      out.par = t_par;
      out.scal = t_scal;
    end
    function out = tperp(obj,species)
      tfac = obj.t_fac(species);
      out = tfac.perp;       
    end
    function out = tpar(obj,species)
      tfac = obj.t_fac(species);
      out = tfac.par;       
    end
    function out = txx(obj,species)
      n = obj.n(species);
      p = obj.pxx(species);
      out = p./n;      
    end
    function out = tyy(obj,species)
      n = obj.n(species);
      p = obj.pyy(species);
      out = p./n;      
    end
    function out = tzz(obj,species)
      n = obj.n(species);
      p = obj.pzz(species);
      out = p./n;      
    end
    function out = tepar(obj)
      species = find(obj.get_charge == -1); % negatively charge particles are electrons 
      tfac = obj.t_fac(species);
      out = tfac.par;       
    end
    function out = teperp(obj)
      species = find(obj.get_charge == -1); % negatively charge particles are electrons 
      tfac = obj.t_fac(species);
      out = tfac.perp;       
    end
    function out = tipar(obj)
      species = find(obj.get_charge == 1); % negatively charge particles are electrons 
      tfac = obj.t_fac(species);
      out = tfac.par;       
    end
    function out = tiperp(obj)
      species = find(obj.get_charge == 1); % negatively charge particles are electrons 
      tfac = obj.t_fac(species);
      out = tfac.perp;       
    end
    
    % Sets of moments, not implemented for Smilei
    function [vxx,vxy,vxz,vyy,vyz,vzz] = vv(obj,iSpecies_orig)
      % [vxx,vxy,vxz,vyy,vyz,vzz] = vv(obj,iSpecies_orig)
      %nargout      
      mass_sp = obj.mass; 
      if numel(unique(mass_sp(iSpecies_orig))) > 1
        error('All species must have the same mass.'); 
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
    function [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = njp(obj,iSpecies,varargin)
      % [n,jx,jy,jz,pxx,pxy,pxz,pyy,pyz,pzz] = njp(obj,iSpecies)
      doMean = 0;
      iSpecies_orig = iSpecies;
      nSpecies = numel(iSpecies);
      dfac = obj.get_dfac;
      args = varargin;
      nargs = numel(args);
      if nargs > 0
        doMean = 1;
        dirMean = args{2};
      end
      
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
        
        n = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        jx = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        jy = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        jz = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vxs = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vys = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vzs = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        pxx = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        pyy = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        pzz = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        pxy = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        pxz = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        pyz = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vxx = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vyy = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vzz = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vxy = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vxz = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        vyz = zeros(numel(obj.xi),numel(obj.zi),obj.length);
        
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
      if doMean
        n = mean(n,dirMean);
        jx = mean(jx,dirMean);
        jy = mean(jy,dirMean);
        jz = mean(jz,dirMean);
        pxx = mean(pxx,dirMean);
        pxy = mean(pxy,dirMean);
        pxz = mean(pxz,dirMean);
        pyy = mean(pyy,dirMean);
        pyz = mean(pyz,dirMean);
        pzz = mean(pzz,dirMean);
      end      
    end  
    
    % Get derived quantities
    function out = minusBx(obj)
      % just for testing obj.get_str('Bx.*Bz');
      out  = obj.get_exp('-Bx');
    end
    function out = par(obj,field,species)
      if exist('species','var') % 
        x = obj.([field ,'x'])(species);
        y = obj.([field ,'y'])(species);
        z = obj.([field ,'z'])(species);
      else % E, B, J
        x = obj.([field ,'x']);
        y = obj.([field ,'y']);
        z = obj.([field ,'z']);
      end
      Bx = obj.Bx;
      By = obj.By;
      Bz = obj.Bz;
      Babs = sqrt(Bx.^2 + By.^2 + Bz.^2);
      bx = Bx./Babs;
      by = By./Babs;
      bz = Bz./Babs;
      out = x.*bx + y.*by + z.*bz;
    end
    function out = parx(obj,field,species)
      if exist('species','var')
        x = obj.([field ,'x'])(species);
        y = obj.([field ,'y'])(species);
        z = obj.([field ,'z'])(species);
      else
        x = obj.([field ,'x']);
        y = obj.([field ,'y']);
        z = obj.([field ,'z']);
      end
      Bx = obj.Bx;
      By = obj.By;
      Bz = obj.Bz;
      Babs = sqrt(Bx.^2 + By.^2 + Bz.^2);
      bx = Bx./Babs;
      by = By./Babs;
      bz = Bz./Babs;
      out = x.*bx;
    end
    function out = pary(obj,field,species)
      if exist('species','var')
        x = obj.([field ,'x'])(species);
        y = obj.([field ,'y'])(species);
        z = obj.([field ,'z'])(species);
      else
        x = obj.([field ,'x']);
        y = obj.([field ,'y']);
        z = obj.([field ,'z']);
      end
      Bx = obj.Bx;
      By = obj.By;
      Bz = obj.Bz;
      Babs = sqrt(Bx.^2 + By.^2 + Bz.^2);
      bx = Bx./Babs;
      by = By./Babs;
      bz = Bz./Babs;
      out = y.*by;
    end
    function out = parz(obj,field,species)
      if exist('species','var')
        x = obj.([field ,'x'])(species);
        y = obj.([field ,'y'])(species);
        z = obj.([field ,'z'])(species);
      else
        x = obj.([field ,'x']);
        y = obj.([field ,'y']);
        z = obj.([field ,'z']);
      end
      Bx = obj.Bx;
      By = obj.By;
      Bz = obj.Bz;
      Babs = sqrt(Bx.^2 + By.^2 + Bz.^2);
      bx = Bx./Babs;
      by = By./Babs;
      bz = Bz./Babs;
      out = z.*bz;
    end
    function varargout = perp(obj,field,comp,species)
      % PIC.PERP Dont know if this one is correct
      % varargout = perp(obj,field,comp,species)
      if exist('species','var') && not(isempty(species))
        x = obj.([field ,'x'])(species);
        y = obj.([field ,'y'])(species);
        z = obj.([field ,'z'])(species);
      else
        x = obj.([field ,'x']);
        y = obj.([field ,'y']);
        z = obj.([field ,'z']);
      end
      Bx = obj.Bx;
      By = obj.By;
      Bz = obj.Bz;
      Babs = sqrt(Bx.^2 + By.^2 + Bz.^2);
      bx = Bx./Babs;
      by = By./Babs;
      bz = Bz./Babs;
      parx = x.*bx; 
      pary = y.*by;
      parz = z.*bz;
      varargout = {};
      for icomp = 1:numel(comp)
        switch comp(icomp)
          case 'x'
            var = x - parx;
          case 'y'
            var = y - pary;
          case 'z'
            var = z - parz;
        end
        varargout{end+1} = var;
      end
    end
    function out = Epar(obj)
      out = obj.par('E');
    end
    function out = Eparx(obj)
      out = obj.parx('E');
    end
    function out = Epary(obj)
      out = obj.pary('E');
    end
    function out = Eparz(obj)
      out = obj.parz('E');
    end
    function out = Experp(obj)
      out = obj.par('E');
    end
    function out = vpar(obj,species)
      out = obj.par('v',species);
    end 
    function out = vepar(obj)
      iSpecies = find(obj.get_charge == -1); % negatively charge particles are electrons
      out = obj.par('v',iSpecies);
    end    
    function out = vipar(obj)
      iSpecies = find(obj.get_charge == 1); % negatively charge particles are electrons
      out = obj.par('v',iSpecies);
    end
    function out = jepar(obj)
      iSpecies = find(obj.get_charge == -1); % negatively charge particles are electrons
      out = obj.par('j',iSpecies);
    end    
    function out = jipar(obj)
      iSpecies = find(obj.get_charge == 1); % negatively charge particles are electrons
      out = obj.par('j',iSpecies);
    end
    function out = dvdt(obj,comp,species)
      % How to deal with dt when nt == 1, i added a parent property to the 
      % PIC class. In that way one can have access to original properties.
      doFirst = 0;
      doLast = 0;
      
      child_it = obj.it;
      child_ix = obj.ix;
      child_iz = obj.iz;
      parent_it = obj.parent.it;
      child_it_required = child_it;
      if not(child_it(1)==1) % add previous timestep
        child_it_required = [child_it_required(1)-1 child_it_required];
      else % otherwise do 2 point derivative for first point
        doFirst = 1;      
      end
      if not(child_it(end)==parent_it(end)) % add following timestep
        child_it_required = [child_it_required child_it_required(end)+1];
      else % otherwise do 2 point derivative for last point
        doLast = 1;
      end
      
      new_obj = obj.parent.subset('t',child_it_required).subset('x',child_ix).subset('z',child_iz);
      
      switch comp
        case 'x', v = new_obj.vx(species);
        case 'y', v = new_obj.vy(species);
        case 'z', v = new_obj.vz(species);
      end
      t = new_obj.twci;     
      dvdt = zeros(size(v));
      nt = new_obj.nt;
      for it = 2:(nt-1)
        % dt might be variable
        dvdt_minus = (v(:,:,it)-v(:,:,it-1))/(t(it)-t(it-1));
        dvdt_plus = (v(:,:,it+1)-v(:,:,it))/(t(it+1)-t(it));
        dvdt(:,:,it) = (dvdt_minus+dvdt_plus)/2;        
      end
      dvdt(:,:,1) = (v(:,:,2)-v(:,:,1))/(t(2)-t(1));
      dvdt(:,:,nt) = (v(:,:,nt)-v(:,:,nt-1))/(t(nt)-t(nt-1));        
      out = dvdt;
      if not(doFirst)
        out = out(:,:,2:end);
      end
      if not(doLast)
        out = out(:,:,1:(end-1));
      end
      
    end
    function out = dvxdt(obj,species)
      out = obj.dvdt('x',species);
    end
    function out = dvydt(obj,species)
      out = obj.dvdt('y',species);
    end
    function out = dvzdt(obj,species)
      out = obj.dvdt('z',species);
    end
    function out = dvixdt(obj)
      species = find(obj.get_charge == 1);
      out = obj.dvdt('x',species);
    end
    function out = dviydt(obj)
      species = find(obj.get_charge == 1);
      out = obj.dvdt('y',species);
    end
    function out = dvizdt(obj)
      species = find(obj.get_charge == 1);
      out = obj.dvdt('z',species);
    end
    function out = dvexdt(obj)
      species = find(obj.get_charge == -1);
      out = obj.dvdt('x',species);
    end
    function out = dveydt(obj)
      species = find(obj.get_charge == -1);
      out = obj.dvdt('y',species);
    end
    function out = dvezdt(obj)
      species = find(obj.get_charge == -1);
      out = obj.dvdt('z',species);
    end
    function out = vdv(obj,comp,species)
      % 
      % vdivv_x = vx*dvxdx + vy*dvydx + vz*dvzdx
      %         = vx*dvxdx + 0        + vz*dvzdx
                  
      dx = obj.xi(2) - obj.xi(1);
      dz = obj.zi(2) - obj.zi(1);
      vx = obj.vx(species);
      vy = obj.vy(species);
      vz = obj.vz(species);
      
      switch comp
        case 'x', v = vx;
        case 'y', v = vy;
        case 'z', v = vz;
      end
      

      [dzV,dxV] = gradient(v,dz,dx);
      dyV = 0;
                   
      vdivv = (dxV + dyV + dzV).*v;
      out = vdivv;
      
    end
    function out = vdvx(obj,species)
      out = obj.vdv('x',species);
    end
    function out = vdvy(obj,species)
      out = obj.vdv('y',species);
    end
    function out = vdvz(obj,species)
      out = obj.vdv('z',species);
    end
    function out = vdvex(obj)
      species = find(obj.get_charge == -1);
      out = obj.vdv('x',species);
    end
    function out = vdvey(obj)
      species = find(obj.get_charge == -1);
      out = obj.vdv('y',species);
    end
    function out = vdvez(obj)
      species = find(obj.get_charge == -1);
      out = obj.vdv('z',species);
    end
    function out = vdvix(obj)
      species = find(obj.get_charge == 1);
      out = obj.vdv('x',species);
    end
    function out = vdviy(obj)
      species = find(obj.get_charge == 1);
      out = obj.vdv('y',species);
    end
    function out = vdviz(obj)
      species = find(obj.get_charge == 1);
      out = obj.vdv('z',species);
    end
    function out = JxBx(obj)
      % vxB_x(vyBz - vzBy)
      % out = vxBx(obj,species)
      By = obj.By;
      Bz = obj.Bz;
      Jy = obj.Jy;
      Jz = obj.Jz;
      out = (Jy.*Bz-Jz.*By);
    end
    function out = JxBy(obj)
      % vxB_y(vzBx - vxBz)
      % out = vxBx(obj,species)
      Bx = obj.Bx;
      Bz = obj.Bz;
      Jx = obj.Jx;
      Jz = obj.Jz;   
      out = (Jz.*Bx-Jx.*Bz);
    end
    function out = JxBz(obj)
      % vxB_z(vxBy - vyBz)
      % out = vxBx(obj,species)
      Bx = obj.Bx;
      By = obj.By;
      Jx = obj.Jx;
      Jy = obj.Jy;
      out = (Jx.*By-Jy.*Bx);
    end
    
    function out = vxBx(obj,species)
      % vxB_x(vyBz - vzBy)
      % out = vxBx(obj,species)
      By = obj.By;
      Bz = obj.Bz;
      vy = obj.vy(species);
      vz = obj.vz(species);      
      out = (vy.*Bz-vz.*By);
    end
    function out = vxBy(obj,species)
      % vxB_y(vzBx - vxBz)
      % out = vxBx(obj,species)
      Bx = obj.Bx;
      Bz = obj.Bz;
      vx = obj.vx(species);
      vz = obj.vz(species);      
      out = (vz.*Bx-vx.*Bz);
    end
    function out = vxBz(obj,species)
      % vxB_z(vxBy - vyBz)
      % out = vxBx(obj,species)
      Bx = obj.Bx;
      By = obj.By;
      vx = obj.vx(species);
      vy = obj.vy(species);
      out = (vx.*By-vy.*Bx);
    end
    function out = vixBx(obj)
      iSpecies = find(obj.get_charge == 1);
      out = obj.vxBx(iSpecies);
    end
    function out = vixBy(obj)
      iSpecies = find(obj.get_charge == 1);
      out = obj.vxBy(iSpecies);
    end
    function out = vixBz(obj)
      iSpecies = find(obj.get_charge == 1);
      out = obj.vxBz(iSpecies);
    end
    function out = vexBx(obj)
      iSpecies = find(obj.get_charge == -1);
      out = obj.vxBx(iSpecies);
    end
    function out = vexBy(obj)
      iSpecies = find(obj.get_charge == -1);
      out = obj.vxBy(iSpecies);
    end
    function out = vexBz(obj)
      iSpecies = find(obj.get_charge == -1);
      out = obj.vxBz(iSpecies);
    end
    function out = vExBx(obj)
      % (EyBz - EzBy)/|B|^2
      Bx = obj.Bx;
      By = obj.By;
      Bz = obj.Bz;
      
      Ey = obj.Ey;
      Ez = obj.Ez;
      
      out = (Ey.*Bz-Ez.*By)./(Bx.^2+By.^2+Bz.^2);
    end
    function out = vExBy(obj)
      % (EyBz - EzBy)/|B|^2
      Bx = obj.Bx;
      By = obj.By;
      Bz = obj.Bz;
      
      Ex = obj.Ex;
      Ez = obj.Ez;
      
      out = (Ez.*Bx-Ex.*Bz)./(Bx.^2+By.^2+Bz.^2);
    end
    function out = vExBz(obj)
      % (EyBz - EzBy)/|B|^2
      Bx = obj.Bx;
      By = obj.By;
      Bz = obj.Bz;
      
      Ex = obj.Ex;
      Ey = obj.Ey;
      
      out = (Ex.*By-Ey.*Bx)./(Bx.^2+By.^2+Bz.^2);
    end
    function out = vt(obj,species,mode) % thermal speed
      % PIC.VT Thermal speed (2T/m)^0.5
      % vt = vt(obj,species,mode)
      % vt = obj.vt(species,mode)
      % To implement: full tensor
      if not(exist('mode','var'))
        mode = 'total';        
      end
      switch mode
        case {'xx','yy','zz'}
          p = obj.p12(species,mode);
          n = obj.n(species);
          t = p./n;
        case {'total','abs','scalar'}
          t = obj.t(species);
        case {'par','perp'}
          tfac = obj.t_fac(species);
          switch mode
            case 'par'
              t = tfac.par;
            case 'perp'
              t = tfac.perp;
          end
      end
      t_neg = find(t<0);
      warning('t has %g negative values, putting them to zero.',numel(t_neg))
      t(t_neg) = 0;
      mass = obj.mass(species(1))/obj.mass(1);
      vt = sqrt(2*t/mass);
      out = vt;
    end
    function out = vtpar(obj,species)
      out = vt(obj,species,'par');
    end
    function out = vtperp(obj,species)
      out = vt(obj,species,'perp');
    end
    function out = vtxx(obj,species)
      out = vt(obj,species,'xx');
    end
    function out = vtyy(obj,species)
      out = vt(obj,species,'yy');
    end
    function out = vtzz(obj,species)
      out = vt(obj,species,'zz');
    end
    function out = wc(obj,species) % plasma frequency
      % PIC.LDE Cyclotron frequency normalized to wp0, should be normalized to
      %   wc0 I suppose since time is normalized to wc0^-1. 
      %   wpe/wce = 2, wci/wce = me/mi, wpi/wpe = sqrt(me/mi)
      %   (wpe/wce)*(wce/wci)*(wpi/wpe) = wpe/wce*mi/me*sqrt(me/mi) =
      %   wpe/wce*sqrt(mi/me) = wpi/wci
      %
      %   out = wc(obj,species);
      %   out = obj.wp(species);
      % 
      Babs = obj.Babs;   
      out = Babs/(obj.mass(species(1))/obj.mass(1));
    end
    function out = rc(obj,species,comp)
      %
      % wc = eB/m
      if not(exist('mode','var'))
        mode = 'total';        
      end
      vt = obj.vt(species,mode);
      wc = obj.wc(species);
      out = vt./wc; % because we have included a 2 inside vt = sqrt(2T/m)
    end
    function out = wp(obj,species) % plasma frequency
      % PIC.LDE Plasma frequency normalized to wp0, should be normalized to
      %   wc0 I suppose since time is normalized to wc0^-1. 
      %   wpe/wce = 2, wci/wce = me/mi, wpi/wpe = sqrt(me/mi)
      %   (wpe/wce)*(wce/wci)*(wpi/wpe) = wpe/wce*mi/me*sqrt(me/mi) =
      %   wpe/wce*sqrt(mi/me) = wpi/wci
      %
      %   out = wp(obj,species);
      %   out = obj.wp(species);
      % 
      n = obj.n(species);
      % wp = sqrt(ne^2/eps0m)
      % should be normalized to wci
      % (wp/wci)^2 = (n*e^2/eps0*m)/(e^2*B0^2/m(1)^2)
      % wp/wp0 = sqrt(n/m)/sqrt(n0/m0) = sqrt(n/n0*m0/m)
      wpwp0 = sqrt(n/obj.mass(species(1)))/sqrt(1/obj.mass(1));
      wpiwci = obj.wpewce*sqrt(obj.mime);
      out = wpwp0*wpiwci;
    end
    function out = lde(obj,species,mode)
      % PIC.LDE Debye length
      %   out = lde(obj,species,mode)
      if not(exist('mode','var'))
        mode = 'total';        
      end
      vt = obj.vt(species,mode);
      wp = obj.wp(species);
      out = vt./wp/sqrt(2); % because we have included a 2 inside vt = sqrt(2T/m)
    end
    function out = ldexx(obj,species)
      out = lde(obj,species,'xx');
    end
    function out = ldeyy(obj,species)
      out = lde(obj,species,'yy');
    end
    function out = ldezz(obj,species)
      out = lde(obj,species,'zz');
    end    
    function out = magmom(obj,species)
      vtperp = obj.vt(species,'perp');
      Bx = obj.Bx;
      By = obj.By;
      Bz = obj.Bz;
      B = sqrt(Bx.^2 + By.^2 + Bz.^2);
      mass = obj.mass(species(1))/obj.mass(1);
      out = 0.5*mass.*vtperp.^2./B;
    end
    function out = gradp(obj,species)
      pxx = obj.pxx(species);
      pyy = obj.pyy(species);
      pzz = obj.pzz(species);
      p = (pxx+pyy+pzz)/3;
      T = p;
            
      dx = obj.xi(2)-obj.xi(1);
      dy = Inf;
      dz = obj.zi(2)-obj.zi(1);

      diff_order = 1; % using diff_order = 2 seems to not work as well
      if diff_order == 1
        dxT = [diff(T(1:2,:),1,1); T(2:end,:)-T(1:(end-1),:)];
        dyT = T*0;
        %dzT = [diff(T(:,1:2),1,2), diff(T,diff_order,2)];
        dzT = [diff(T(:,1:2),1,2), T(:,2:end)-T(:,1:(end-1))];
      else diff_order == 2
        dxT = [2*diff(T(1:2,:),1,1); diff(T,diff_order,1); 2*diff(T((end-1):end,:),1,1)];
        dyT = T*0;
        dzT = [2*diff(T(:,1:2),1,2), diff(T,diff_order,2), 2*diff(T(:,(end-1):end),1,2)];
      end

      div_x = dxT/dx/diff_order;
      div_y = dyT/dy/diff_order;
      div_z = dzT/dz/diff_order;

      out.x = div_x;
      out.y = div_y;
      out.z = div_z;
      out.abs = sqrt(div_x.^2 + div_y.^2 + div_z.^2);
    end
    function out = gradpx(obj,species)
      gradp = obj.gradp(species);
      out = gradp.x;
    end
    function out = gradpz(obj,species)
      gradp = obj.gradp(species);
      out = gradp.z;
    end    
    function varargout = divp(obj,species,varargin)
      returnComponents = 0;
      diff_order = 1;
      np_smooth = 7;
      doSmooth = 0;
      nargs = numel(varargin);
      if nargs > 1
        if strcmp(varargin{1},'comp')
          comp = varargin{2};
        end
      end

      nx = obj.nx;
      nz = obj.nz;

      dx = obj.xi(2)-obj.xi(1);
      dy = Inf;
      dz = obj.zi(2)-obj.zi(1);
      varargout = cell(0);
      for icomp = 1:numel(comp)
        switch comp{icomp}
          case {'x',1}
            pxx = obj.pxx(species);
            pxy = obj.pxy(species);
            pxz = obj.pxz(species);
            dxTxx = [1*diff(pxx(1:2,:),1,1); diff(pxx,diff_order,1)]/dx/diff_order;
            dyTyx = 0;
            dzTzx = [1*diff(pxz(:,1:2),1,2), diff(pxz,diff_order,2)]/dz/diff_order;
            div_x = dxTxx + dyTyx + dzTzx;
            varargout{end+1} = div_x;
          case {'y',2}
            % pyy = obj.pxx(species); % not needed because dy = inf;
            pxy = obj.pxy(species);
            pyz = obj.pyz(species);
            if 0
              dxTxy = [1*diff(pxy(1:2,:),1,1); diff(pxy,diff_order,1)]/dx/diff_order;
              dyTyy = 0;
              dzTzy = [1*diff(pyz(:,1:2),1,2), diff(pyz,diff_order,2)]/dz/diff_order;              
            else
%               Note: The first output FX is always the gradient along the 2nd
%               dimension of F, going across columns.  The second output FY is always
%               the gradient along the 1st dimension of F, going across rows.  For the
%               third output FZ and the outputs that follow, the Nth output is the
%               gradient along the Nth dimension of F.
%               [dxTxy,dzTxy] = gradient(pxy,dx,dz);
%               dyTyy = 0;
%               [dxTzy,dzTzy] = gradient(pyz,dx,dz);
              [dzTxy,dxTxy] = gradient(pxy,dz,dx);
              dyTyy = 0;
              [dzTzy,dxTzy] = gradient(pyz,dz,dx);
            end
            div_y = (dxTxy + dyTyy + dzTzy);
            varargout{end+1} = div_y;
          case {'z',3}
            pzz = obj.pzz(species);
            pxz = obj.pxz(species);
            pyz = obj.pyz(species);
            dxTxz = [1*diff(pxz(1:2,:),1,1); diff(pxz,diff_order,1)]/dx/diff_order;
            dyTyz = 0;
            dzTzz = [1*diff(pzz(:,1:2),1,2), diff(pzz,diff_order,2)]/dz/diff_order;      
            div_z = (dxTxz + dyTyz + dzTzz);
            varargout{end+1} = div_z;
        end
      end
      if doSmooth
        for ivar = 1:numel(varargout)
          varargout{ivar} = smooth2(varargout{ivar},np_smooth,np_smooth);
        end
      end

      % dxTxx = zeros(nx,nz);
      % dxTxy = zeros(nx,nz);
      % dxTxz = zeros(nx,nz);
      % dyTyx = zeros(nx,nz);
      % dyTyy = zeros(nx,nz);
      % dyTyz = zeros(nx,nz);
      % dzTzx = zeros(nx,nz);
      % dzTzy = zeros(nx,nz);
      % dzTzz = zeros(nx,nz);

      %out.x = div_x;
      %out.y = div_y;
      %out.z = div_z;

      if 0%returnComponents % is this the correct placements ???
        out.x_xx = dxTxx;
        out.x_yy = dyTyx;
        out.x_zz = dzTzx;
        out.y_xx = dxTxy;
        out.y_yy = dyTyy;
        out.y_zz = dzTzy;
        out.z_xx = dxTxz;
        out.z_yy = dyTyz;
        out.z_zz = dzTzz;
      end
    end
    function out = divpx(obj,species)
      out = obj.divp(species,'comp',{'x'});
    end
    function out = divpy(obj,species)
      out = obj.divp(species,'comp',{'y'});
    end
    function out = divpz(obj,species)
      out = obj.divp(species,'comp',{'z'});
    end
    function out = divpix(obj)
      species = find(obj.get_charge == 1);
      out = obj.divp(species,'comp',{'x'});
    end
    function out = divpiy(obj)
      species = find(obj.get_charge == 1);
      out = obj.divp(species,'comp',{'y'});
    end
    function out = divpiz(obj)
      species = find(obj.get_charge == 1);
      out = obj.divp(species,'comp',{'z'});
    end
    function out = divpex(obj)
      species = find(obj.get_charge == -1);
      out = obj.divp(species,'comp',{'x'});
    end
    function out = divpey(obj)
      species = find(obj.get_charge == -1);
      out = obj.divp(species,'comp',{'y'});
    end
    function out = divpez(obj)
      species = find(obj.get_charge == -1);
      out = obj.divp(species,'comp',{'z'});
    end
    % Stored, not implemented for Smilei
    function out = UB(obj,varargin)
      % Magnetic energy density 0.5*(Bx^2 + By^2 + Bz^2) summed up
      out = obj.get_timeline_attributes('UB');
      if numel(varargin)>0 && isnumeric(varargin{1})
        inds = varargin{1};
        out = out(inds);
      end
      %out = h5read(obj.file,'/scalar_timeseries/U/B');
      %out = out(obj.indices_);
    end
    function out = dUB(obj)
      % Magnetic energy density 0.5*(Bx^2 + By^2 + Bz^2) summed up 
      out = obj.UB';
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
    function out = Ute(obj)
      % Integrated electron thermal energy density
      out = obj.get_timeline_attributes('Ute');
    end
    function out = Uti(obj)
      % Integrated ion thermal energy density
      out = obj.get_timeline_attributes('Uti');
    end
    function out = Uke(obj)
      % Integrated electron kinetic energy density
      out = obj.get_timeline_attributes('Uke');
    end
    function out = Uki(obj)
      % Integrated electron kinetic energy density
      out = obj.get_timeline_attributes('Uki');
    end
    function out = RE(obj)
      % Reconnection rate from out-of-plane electric field Ey at X line
      try
        out = obj.get_timeline_attributes('RE');
      catch
        try
          out = h5read(obj.file,'/scalar_timeseries/R/Ey');
          out = out(obj.indices_);
        end
      end
    end
    function out = RA(obj)
      % Reconnection rate from vector potential dA/dt at X line
      % out = h5read(obj.file,'/scalar_timeseries/R/A');      
      out = obj.get_timeline_attributes('RA');        
      if all(isnan(out))
        disp('''out = obj.get_timeline_attributes(''RA'');'' did not work.')
        disp('Trying: ''out = obj.get_timeline_attributes(''Axline'');''')      
        Aval = obj.get_timeline_attributes('Axline');
        dA = diff(Aval);
        dt = diff(obj.twci);
        times = obj.twci;
        dAdt_ = dA./dt;
        dAdt = interp1(times(1:end-1)+0.5*dt,dAdt_,times);
        out  = -dAdt;
      end
      out = out(obj.indices_);
    end
    function out = Axline(obj)
      % A value at X line position
      A_xline = obj.get_timeline_attributes('Axline');
      out = A_xline(obj.indices,:);
    end
    function out = xline_position(obj)
      % X line position
      xz_xline = obj.get_timeline_attributes('xline_position');
      out = xz_xline(obj.indices,:);
    end
    function out = x_xline(obj)
      % X line position
      xz_xline = obj.get_timeline_attributes('xline_position');
      out = xz_xline(obj.indices_,1);
    end
    function out = z_xline(obj)
      % X line position
      xz_xline = obj.get_timeline_attributes('xline_position');
      out = xz_xline(obj.indices_,2);
    end
    % Calculated each time, not implemented for Smilei
    function out = curlb(obj)
      % PIC.CURLB Rotation of B.
      %   Jx = dBy/dz - dBz/dy
      %   Jy = dBz/dx - dBx/dz
      %   Jz = dBx/dy - dBy/dx
      %
      %   h = setup_subplots(3,3);
      %
      %
            
      bx = obj.Bx;
      by = obj.By;
      bz = obj.Bz;
      
      dx = obj.xi(2)-obj.xi(1);
      dz = obj.zi(2)-obj.zi(1);
      dy = Inf;
      
      
      % 2D mesh of grid points, for interpolating
      [XI_EY,ZI_EY] = meshgrid(obj.xivar.Ey,obj.zivar.Ey);
      [XI_BY,ZI_BY] = meshgrid(obj.xivar.By,obj.zivar.By);
      [XI_BX,ZI_BX] = meshgrid(obj.xivar.Bx,obj.zivar.Bx);
      [XI_BZ,ZI_BZ] = meshgrid(obj.xivar.Bz,obj.zivar.Bz);
      
      % Since we are not always at the edge of the box, I just treat the
      % edges as duplicates of the adjacent value. 
      
      % -- adjusted grid
      dxbz = zeros(obj.nx,obj.nz);
      dxbz(2:end,:) = diff(bz,1,1)/dx;
      dxbz(1,:) = dxbz(2,1);
      dybz = 0;
      
      dzbx = zeros(obj.nx,obj.nz);
      dzbx(:,1) = (bx(:,2)-bx(:,1))/(dz); % just the same as the pount inside
      dzbx(:,2:end) = diff(bx,1,2)/dz;
      dybx = 0;
      
      dxby_ = diff(by,1,1)/dx; % ends up on Bx/Ez grid one step to the right
      dxby = interp2(XI_BX(:,2:end),ZI_BY(:,2:end),dxby_',XI_EY,ZI_EY)';
      
      dzby_ = diff(by,1,2)/dz; % ends up on Bz/Ex grid one step to the top
      dzby = interp2(XI_BY(2:end,:),ZI_BY(2:end,:),dzby_',XI_EY,ZI_EY)';
      
      %bx_at_ey = interp2(XI_BX,ZI_BX,bx',XI_EY,ZI_EY)';
      %by_at_ey = interp2(XI_BY,ZI_BY,by',XI_EY,ZI_EY)';
      %bz_at_ey = interp2(XI_BZ,ZI_BZ,bz',XI_EY,ZI_EY)';
      
      % Jx = dBy/dz - dBz/dy
      % Jy = dBz/dx - dBx/dz
      % Jz = dBx/dy - dBy/dx
      bcurl_x = dybz - dzby;
      bcurl_y = dzbx - dxbz;
      bcurl_z = dxby - dybx;
      
      out.x = bcurl_x;
      out.y = bcurl_y;
      out.z = bcurl_z;      
    end
    function out = magnetic_curvature(obj)
      % PIC.MAGNETIC_CURVATURE Calculate magnetic field curvature.
      %   bcurv = dot(b,nabla) b
      %   bcurv_x = (bxdx + bydy + bzdz) bx = bx*(dxbx) + by*(dybx) + bz*(dzbx)
      %   bcurv_y = (bxdx + bydy + bzdz) by = bx*(dxby) + by*(dyby) + bz*(dzby)
      %   bcurv_z = (bxdx + bydy + bzdz) bz = bx*(dxbz) + by*(dybz) + bz*(dzbz)
      %   where b is the magnetic field unit vector b = B/|B|
      %
      %   KB = PIC.MAGNETIC_CURVATURE;
      %     KB is centered on "basic grid": pic.xi/pic.zi, pic.xe/pic.ze.
      doParent = 0;
      % Interpolating (I think) makes some end points nan (z end points i think).
      % So use the parent to in
      if obj.iz(1) > 2 && obj.iz(end) < obj.parent.nz-1 && obj.ix(1) > 2 && obj.ix(end) < obj.parent.nx-1
        doParent = 1;
      end
      if doParent
        obj_orig = obj;
        obj = obj_orig.parent;
        obj = obj.twpelim(obj_orig.twpe);
        old_inds_x = obj_orig.ix;
        old_inds_z = obj_orig.iz;
        new_inds_x = obj_orig.ix; new_inds_x = [new_inds_x(1)+[-2 -1] new_inds_x new_inds_x(end)+[1 2]];
        new_inds_z = obj_orig.iz; new_inds_z = [new_inds_z(1)+[-2 -1] new_inds_z new_inds_z(end)+[1 2]];
        new_inds_x = intersect(obj.ix,new_inds_x); % remove if any are outside max
        new_inds_z = intersect(obj.iz,new_inds_z); % remove if any are outside max      
        obj = obj.subset('x',new_inds_x).subset('z',new_inds_z);
        diff_ind_x = setdiff(new_inds_x,old_inds_x);
        diff_ind_z = setdiff(new_inds_z,old_inds_z);
      end
      
      % Load B, consider loading one extra point
      Bx = obj.Bx;
      By = obj.By;
      Bz = obj.Bz;
      
      % Calculate b = B/|B|
      Babs = sqrt(Bx.^2 + By.^2 + Bz.^2);
      bx = Bx./Babs;
      by = By./Babs;
      bz = Bz./Babs;

      % Grid spacing, dx, dy, dz
      dx = obj.xi(2)-obj.xi(1);
      dz = obj.zi(2)-obj.zi(1);
      dy = Inf;
      
      % All the magnetic field components are offset from the central grid
      % pic.xi, pic.zi, see pic.x.i_Bx, pic.x.i_By, pic.x.i_Bz:
      % Bz is half a grid space off to the right,
      % Bx is half a grid space off to the top,
      % By is half a grid space off to the top and right.
      % This will make dxBz, dxBy end up on the central x grid, and
      % dzBx, dzBy end up on the central z grid. For dxBx, and dzBz,
      % however, the grid needs to be shifted/interpolated first.
      
      % Boundary conditions are special cases, because we are not always
      % loading the entire box... Just make a duplicate in the empty
      % cell... But choose the cell based on grid. No, just do interpolate
      % instead.
      
%       if 1
      % 2D mesh of grid points, for interpolating
      [XI_EY,ZI_EY] = meshgrid(obj.xivar.Ey,obj.zivar.Ey);
      [XI_BY,ZI_BY] = meshgrid(obj.xivar.By,obj.zivar.By);
      [XI_BX,ZI_BX] = meshgrid(obj.xivar.Bx,obj.zivar.Bx);
      [XI_BZ,ZI_BZ] = meshgrid(obj.xivar.Bz,obj.zivar.Bz);
      
      % -- new
      dxbz = zeros(obj.nx,obj.nz);
      dxbz(2:end,:) = diff(bz,1,1)/dx;
      dxbz(1,:) = dxbz(2,1);
      dybz = 0;
      dzbz_ = diff(bz,1,2); % this ends up on By part of grid, needs to be 
                            % moved half grid to right or left and half up
                            % or down.
      % interpolating from By grid points to Ey grid point
      dzbz = interp2(XI_BY(1:end-1,:),ZI_BY(1:end-1,:),dzbz_',XI_EY,ZI_EY)';
      
      dzbx = zeros(obj.nx,obj.nz);
      dzbx(:,1) = (bx(:,2)-bx(:,1))/(dz); % just the same as the point inside
      dzbx(:,2:end) = diff(bx,1,2)/dz;
      dybx = 0;
      dxbx_ = diff(bx,1,2); % this ends up on By part of grid
      % interpolating from By grid points to Ey grid point
      dxbx = interp2(XI_BY(1:end-1,:),ZI_BY(1:end-1,:),dxbx_',XI_EY,ZI_EY)';
      
      dxby_ = diff(by,1,1)/dx; % ends up on Bx/Ez grid one step to the right
      dxby = interp2(XI_BX(:,2:end),ZI_BY(:,2:end),dxby_',XI_EY,ZI_EY)';
      
      dzby_ = diff(by,1,2)/dz; % ends up on Bz/Ex grid one step to the top
      dzby = interp2(XI_BY(2:end,:),ZI_BY(2:end,:),dzby_',XI_EY,ZI_EY)';
      dyby = 0;
      
      bx_at_ey = interp2(XI_BX,ZI_BX,bx',XI_EY,ZI_EY)';
      by_at_ey = interp2(XI_BY,ZI_BY,by',XI_EY,ZI_EY)';
      bz_at_ey = interp2(XI_BZ,ZI_BZ,bz',XI_EY,ZI_EY)';
      
      % quick fix for nan endpoints due to interp2, assign closest value
      bx_at_ey(:,1) = bx_at_ey(:,2); bx_at_ey(:,end) = bx_at_ey(:,end-1);
      by_at_ey(:,1) = by_at_ey(:,2); by_at_ey(:,end) = by_at_ey(:,end-1);      
      bz_at_ey(:,1) = bz_at_ey(:,2); bz_at_ey(:,end) = bz_at_ey(:,end-1);
      bx_at_ey(1,:) = bx_at_ey(2,:); bx_at_ey(end,:) = bx_at_ey(end-1,:);
      by_at_ey(1,:) = by_at_ey(2,:); by_at_ey(end,:) = by_at_ey(end-1,:);
      bz_at_ey(1,:) = bz_at_ey(2,:); bz_at_ey(end,:) = bz_at_ey(end-1,:);
      dxbx(:,1) = dxbx(:,2); dxbx(:,end) = dxbx(:,end-1);
      dxbx(1,:) = dxbx(2,:); dxbx(end,:) = dxbx(end-1,:);
      dzbz(:,1) = dzbz(:,2); dzbz(:,end) = dzbz(:,end-1);
      dzbz(1,:) = dzbz(2,:); dzbz(end,:) = dzbz(end-1,:);
      dxby(:,1) = dxby(:,2); dxby(:,end) = dxby(:,end-1);
      dxby(1,:) = dxby(2,:);
      dzby(:,1) = dzby(:,2); dzby(:,end) = dzby(:,end-1);
      dzby(:,2) = dzby(:,3);
      dzby(:,1) = dzby(:,2); dzby(:,end) = dzby(:,end-1);
      dzby(1,:) = dzby(2,:); dzby(end,:) = dzby(end-1,:);
      
      bcurv_x = bx_at_ey.*dxbx + by_at_ey.*dybx + bz_at_ey.*dzbx;
      bcurv_y = bx_at_ey.*dxby + by_at_ey.*dyby + bz_at_ey.*dzby;
      bcurv_z = bx_at_ey.*dxbz + by_at_ey.*dybz + bz_at_ey.*dzbz;
      
      1;
%       else % -- old
%       dxbx = diff(bx,1,1)/dx; dxbx(end+1,:) = dxbx(end,:);
%       dybx = 0;
%       dzbx = diff(bx,1,2); dzbx(:,end+1) = dzbx(:,end);
%       dxby = diff(by,1,1); dxby(end+1,:) = dxby(end,:);
%       dyby = 0;
%       dzby = diff(by,1,2); dzby(:,end+1) = dzby(:,end);
%       dxbz = diff(bz,1,1); dxbz(end+1,:) = dxbz(end,:);
%       dybz = 0;
%       dzbz = diff(bz,1,2); dzbz(:,end+1) = dzbz(:,end);
%       bcurv_x = bx.*dxbx/dx + by.*dybx/dy + bz.*dzbx/dz;
%       bcurv_y = bx.*dxby/dx + by.*dyby/dy + bz.*dzby/dz;
%       bcurv_z = bx.*dxbz/dx + by.*dybz/dy + bz.*dzbz/dz;
%       end

      if doParent
        bcurv_x = bcurv_x(3:end-2,3:end-2);
        bcurv_y = bcurv_y(3:end-2,3:end-2);
        bcurv_z = bcurv_z(3:end-2,3:end-2);
      end
      
      
      bcurv.units = 'wpi/c';
      bcurv.x = bcurv_x;
      bcurv.y = bcurv_y;
      bcurv.z = bcurv_z;
      bcurv.abs = sqrt(bcurv_x.^2 + bcurv_y.^2 + bcurv_z.^2);
      out = bcurv;
      
    end
    function out = curvb(obj)
      %
      % See also PIC.MAGNETIC_CURVATURE
      out = obj.magnetic_curvature;
    end
    function out = curvbx(obj)
      %
      % See also PIC.MAGNETIC_CURVATURE
      curvb = obj.curvb;
      out = curvb.x; 
    end
    function out = curvby(obj)
      %
      % See also PIC.MAGNETIC_CURVATURE
      curvb = obj.curvb;
      out = curvb.y; 
    end
    function out = curvbz(obj)
      %
      % See also PIC.MAGNETIC_CURVATURE
      curvb = obj.curvb;
      out = curvb.z; 
    end
    function out = curvbabs(obj)
      %
      % See also PIC.MAGNETIC_CURVATURE
      curvb = obj.curvb;
      out = curvb.abs; 
    end
    function out = curvbrad(obj)
      %
      % See also PIC.MAGNETIC_CURVATURE
      curvb = obj.curvb;
      out = 1./curvb.abs; 
    end
    function out = xline(obj)
      % Not implemented.
      % Assume xline is close to middle of box in z
      % zlim = [-0.2 0.2];
      % Assume xline has not moved too much left and right
      %xlim = [-50 50];
      [saddle_x,saddle_z,saddle_val] = obj.saddle;
      A = obj.A; 
      1;
    end
    function varargout = saddle(obj)
      % Not implemented.
      % Assume xline is close to middle of box in z
      %zlim = [-0.2 0.2];
      % Assume xline has not moved too much left and right
      %xlim = [-50 50];
      %A = obj.A;
      times = obj.twci;
      for it = 1:obj.nt
        A = obj.twcilim(times(it)).A;
        [inds,vals] = saddle(A,'sort');
        inds_x_all(it) = inds(1,1);
        inds_z_all(it) = inds(1,2);
        vals_all(it) = vals(1);
      end
      
      varargout{1} = times;
      varargout{2} = obj.xi(inds_x_all);
      varargout{3} = obj.zi(inds_z_all);
      varargout{4} = vals_all;
    end
    function [x,v,a,B] = xva_df(obj)
      % [xDF,vDF,aDF,BDF] = df04.xva_df;
      % Divided into left and right o center of box
      
      % Find main X line and divide box into left and right of this
      % Never mind, jsut pick 0.
      
      nt = obj.nt;
      tlim = obj.twci([1 end]);
      zlim = [-0.2 0.2];
      xlims = {[obj.xi(1) obj.xi(fix(obj.nx/2))],[obj.xi(fix(obj.nx/2)) obj.xi(end)]};
      
      x = zeros(2,obj.nt);
      v = zeros(2,obj.nt);
      a = zeros(2,obj.nt);
      B = zeros(2,obj.nt);
      multiplier = [-1 1];
      for ixlim = 1:2 
        xlim = xlims{ixlim};        
        pic = obj.xlim(xlim).zlim(zlim).twcilim(tlim);
        dt = reshape(diff(pic.twci),[nt-1,1]);
        t_centered = reshape(pic.twci(1:end-1),[nt-1 1])+dt;
        
        Bz_mean = squeeze(mean(pic.Bz,2)); % mean over z range
        [Bz_peak,ind_Bz_peak] = max(abs(Bz_mean)); % find peak Bz
        xDF = reshape(pic.xi(ind_Bz_peak),[nt,1]); % get x-locations of peak Bz
        
        dxDF = diff(xDF);
        vDF = dxDF./dt;
        vDF_interp = reshape(interp1(t_centered,vDF,pic.twci),[nt,1]); % interpolate to original timeseries
        %tcentered_2 = tocolumn(pic.twci(2:pic.nt-1));
        aDF = reshape(diff(vDF),[nt-2,1])./dt(2:end);
        aDF = [aDF(1); aDF; aDF(end)];     
        x(ixlim,:) = xDF;
        v(ixlim,:) = vDF_interp*multiplier(ixlim);
        a(ixlim,:) = aDF*multiplier(ixlim);
        B(ixlim,:) = Bz_peak;
      end
    end  
    function out = fft(obj,dim)
      
    end
    
    function out = xmesh(obj)
      x = obj.xi;
      out = repmat(x,1,obj.nz);
    end
    function out = zmesh(obj)
      z = obj.zi;
      out = repmat(torow(z),obj.nx,1);
    end
  end
  methods % Get and set properties
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
    function obj = set.ye(obj,value)
      obj.ye_ = value;
    end
    function obj = set.ze(obj,value)
      obj.ze_ = value;
    end
    function obj = set.xi(obj,value)
      obj.xi_ = value;
    end
    function obj = set.yi(obj,value)
      obj.yi_ = value;
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
    function obj = set.iy(obj,value)
      obj.iy_ = value;
    end
    function obj = set.iz(obj,value)
      obj.iz_ = value;
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
    function value = get.ye(obj)
      value = obj.ye_;
    end
    function value = get.ze(obj)
      value = obj.ze_;
    end
    function value = get.xi(obj)
      value = obj.xi_;
    end
    function value = get.yi(obj)
      value = obj.yi_;
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
    function value = get.iy(obj)
      value = obj.iy_;
    end
    function value = get.iz(obj)
      value = obj.iz_;
    end
  end
  methods (Static) % does not require object as input, but still needs to be called as obj.func (?)
    function out = ind_from_lim(var,value,varargin)
      % method is the same for xlim, zlim ilim, i, twpelim, twcilim
      
      % Defaults      
      doExact = 0; % find all exact matches, can be any number
      doBounding = 0; % find a given number of values bounding the given value
      doBoundingOld = 0;
      nBounding = 0;
      doClosest = 0;
      nClosest = 1; % only the closest index, can be one or many, for example to reduce cadence
      
      if numel(value) == 1
        doClosest = 1;        
      end      
      if numel(value) > 2
        doExact = 1;        
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
          case 'exact'
            l = 1;
            doExact = 1;
            args = args(l+1:end);
          otherwise
            warning(sprintf('Input ''%s'' not recognized.',args{1}))
            args = args(l+1:end);
        end        
        if isempty(args), break, end    
      end
      
      % Find indices
      if doBounding
        i1 = find(var<value(1),nBounding,'last');
        i2 = find(var>value(1),nBounding,'first');
        
        % Check so that indices are not outside range
        if i1 < 1
          i1 = 1; 
        end
        if i2 > numel(var) 
          i2 = numel(var); 
        end
        inds = i1:i2;
      elseif doBoundingOld
        i0 = find(abs(var-value(1)) == min(abs(var-value(1))));
        i1 = i0 - nBounding;
        i2 = i0 + nBounding;
        % Check so that indices are not outside range
        if i1 < 1
          i1 = 1; 
        end
        if i2 > numel(var) 
          i2 = numel(var); 
        end
        inds = i1:i2;
      elseif doClosest        
        ii = abs(var-value(1));
        [is, index] = sort(abs(var-value(1)));
        inds = sort(index(1:nClosest)); % why did i just not take the first closest index?     
      elseif doExact
        [~,inds,~] = intersect(var,value);
      else        
        i1 = find(var >= value(1),1,'first'); 
        i2 = find(var <= value(2),1,'last'); 
        inds = i1:i2;
      end
      
      out = inds;
    end
    function out = eom()
    end  
    function vargout = rotate_tens(varargin)
      % new_tens = R*old_tens*R^T

      switch nargin 
        case 2 % constant rotation matrix
          old_xx = varargin{1}.xx;
          old_xy = varargin{1}.xy;
          old_xz = varargin{1}.xz;
          old_yy = varargin{1}.yy;
          old_yz = varargin{1}.yz;
          old_zz = varargin{1}.zz;        
          rx = varargin{2}(1,:); 
          ry = varargin{2}(2,:);
          rz = varargin{2}(3,:);
        case 4 % spatially varying rotation matrix
          old_xx = varargin{1}.xx;
          old_xy = varargin{1}.xy;
          old_xz = varargin{1}.xz;
          old_yy = varargin{1}.yy;
          old_yz = varargin{1}.yz;
          old_zz = varargin{1}.zz;        
          rx = varargin{2};
          ry = varargin{3};
          rz = varargin{4};
        case 9
          old_xx = varargin{1};
          old_xy = varargin{2};
          old_xz = varargin{3};
          old_yy = varargin{4};
          old_yz = varargin{5};
          old_zx = varargin{6};    
          rx = varargin{7};
          ry = varargin{8};
          rz = varargin{9};    
        otherwise
          error('Input not recognized.')    
      end
      %      | rx.x rx.y rx.z |
      % r =  | ry.x ry.y ry.z |
      %      | rz.x rz.y rz.z |
      %
      %      | rx.x ry.x rz.x |
      % rt = | rx.y ry.y rz.y |
      %      | rx.z ry.z rz.z |
      %                       
      % newT = r*T*rt
      %
      %        | rx.x rx.y rx.z |     | T.xx T.xy T.xz |     | rx.x ry.x rz.x |
      %     =  | ry.x ry.y ry.z | dot | T.yx T.yy T.yz | dot | rx.y ry.y rz.y |
      %        | rz.x rz.y rz.z |     | T.zx T.zy T.zz |     | rx.z ry.z rz.z |
      %
      %        | rx.x*T.xx + rx.y*T.yx + rx.z*T.zx |     | rx.x ry.x rz.x |
      %     =  | ry.x*T.xy + ry.y*T.yy + ry.z*T.zy |  dot | rx.y ry.y rz.y |
      %        | rz.x*T.xz + rz.y*T.yz + rz.z*T.zz |      | rx.z ry.z rz.z |
      %
      %        | rTx |'    | rx.x ry.x rz.x |
      %     =  | rTy | dot | rx.y ry.y rz.y |
      %        | rTz |     | rx.z ry.z rz.z |


      %rxt.x = rx.x; rxt.y = ry.x; rxt.z = rz.x;
      %ryt.x = rx.y; ryt.y = ry.y; ryt.z = rz.x;

      R = zeros(3,3,size(old_xx,1),size(old_xx,2));
      R(1,1,:,:) = rx.x;
      R(1,2,:,:) = rx.y;
      R(1,3,:,:) = rx.z;
      R(2,1,:,:) = ry.x;
      R(2,2,:,:) = ry.y;
      R(2,3,:,:) = ry.z;
      R(3,1,:,:) = rz.x;
      R(3,2,:,:) = rz.y;
      R(3,3,:,:) = rz.z;

      T = zeros(3,3,size(old_xx,1),size(old_xx,2));
      T(1,1,:,:) = old_xx;
      T(1,2,:,:) = old_xy;
      T(1,3,:,:) = old_xz;
      T(2,1,:,:) = old_xy;
      T(2,2,:,:) = old_yy;
      T(2,3,:,:) = old_yz;
      T(3,1,:,:) = old_xz;
      T(3,2,:,:) = old_yz;
      T(3,3,:,:) = old_zz;

      newT = T*0;
      % sum over i j
      % nT_mn = R_mi * R_nj * T_ij
      for mm = 1:3
        for nn = mm:3
          for ii = 1:3
            for jj = 1:3
              newT(mm,nn,:,:) = newT(mm,nn,:,:) + R(mm,ii,:,:).*R(nn,jj,:,:).*T(ii,jj,:,:);
            end
          end
        end
      end

      new_xx = squeeze(newT(1,1,:,:));
      new_xy = squeeze(newT(1,2,:,:));
      new_xz = squeeze(newT(1,3,:,:));
      new_yy = squeeze(newT(2,2,:,:));
      new_yz = squeeze(newT(2,3,:,:));
      new_zz = squeeze(newT(3,3,:,:));

      % rTx = old_xx.*rx.x + old_xy.*rx.y + old_xz.*rx.z;
      % rTy = old_xy.*rx.x + old_yy.*rx.y + old_yz.*rx.z;
      % rTz = old_xz.*rx.x + old_yz.*rx.y + old_zz.*rx.z;
      % 
      % new_xx = rTx.*rx.x + rTx.*rx.y + rTx.*rx.z;
      % 
      % new_yy = old_xx.*ry.x + old_xy.*ry.y + old_xz.*ry.z;
      % new_xy = old_xy.*ry.x + old_yy.*ry.y + old_yz.*ry.z;
      % new_xz = old_xz.*rx.x + old_yz.*rx.y + old_zz.*rx.z;
      % 
      % 
      % new_y = old_x.*ry.x + old_y.*ry.y + old_z.*ry.z;
      % new_z = old_x.*rz.x + old_y.*rz.y + old_z.*rz.z;

      switch nargout
        case 1
          new_tens.xx = new_xx;
          new_tens.xy = new_xy;
          new_tens.xz = new_xz;
          new_tens.yy = new_yy;
          new_tens.yz = new_yz;
          new_tens.zz = new_zz;
          vargout(1) = new_tens;
        case 6
          vargout(1) = new_xx;
          vargout(2) = new_xy;
          vargout(3) = new_xz;
          vargout(4) = new_yy;
          vargout(5) = new_yz;
          vargout(6) = new_zz;
      end
    end
    function out = curl_vector(x,z,vx,vy,vz)
    end
    function fout = interpfield3(x,y,z,f,xq,yq,zq)
      % INTERPFIELD3 Linear interpolation of 3D field.
      %   bout = interpfield(x,z,t,f,xq,yq,zq)
      %   x, z, t - grid
      %   f - field to interploate (defined on grid x,z,t)  
      %   xq, zq, tq - point to interpolate to, can be scalar values or vectors
      %
      nq = numel(xq);
      for iq = 1:nq
        % find two closest grid points in x and z
        %   11 ------- 21
        %    |  o      |      o is point (xq,zq)
        %    |         |
        %    |         |
        %   12---------22  

        % First do two 2D interpolations, then one 1D interpolation
        ix1 = find(x<xq(iq),1,'last');
        ix2 = find(x>xq(iq),1,'first');
        iy1 = find(y<yq(iq),1,'last');
        iy2 = find(y>yq(iq),1,'first');
        iz1 = find(z<zq(iq),1,'last');
        iz2 = find(z>zq(iq),1,'first');
        % pick out x and z values for these points
        x1 = x(ix1); % x
        x2 = x(ix2);
        y1 = y(iy1); % z
        y2 = y(iy2); 
        z1 = z(iz1); % t 
        z2 = z(iz2);  
        % pick out field values at the eight corners ()
        f111 = f(ix1,iy1,iz1);
        f112 = f(ix1,iy1,iz2);
        f121 = f(ix1,iy2,iz1);
        f122 = f(ix1,iy2,iz2);
        f211 = f(ix2,iy1,iz1);
        f212 = f(ix2,iy1,iz2);
        f221 = f(ix2,iy2,iz1);
        f222 = f(ix2,iy2,iz2);  
        % distances between point 'o' and edges of box (x1,x2,z1,z2)
        xx = (xq(iq) - x1)/(x2-x1); %
        yy = (yq(iq) - y1)/(y2-y1);  
        zz = (zq(iq) - z1)/(z2-z1);
        % bilinear interpolation
        % you can solve for coefficients by matrix operations, but for now I just
        % used the expression directly: https://en.wikipedia.org/wiki/Bilinear_interpolation
        % Here are instructions for 3D: https://en.wikipedia.org/wiki/Trilinear_interpolation    
        % First do two 2D interpolations in yz-plane, then one 1D interpolation
        % along x (which order you do the planes in doesn't matter)
        % Bilinear
        f1 = f111*(1-xx)*(1-yy) + f121*xx*(1-yy) + f112*(1-xx)*yy + f122*xx*yy;
        f2 = f211*(1-xx)*(1-yy) + f221*xx*(1-yy) + f212*(1-xx)*yy + f222*xx*yy;
        % Linear
        f12(iq) = f1*(1-zz) + f2*zz;

        % debug
        %disp(sprintf('istep = %g, [x1,x1,z1,z1] = [%g,%g,%g,%g],  [b11,b12,b21,b22] = [%g,%g,%g,%g], binterp = %g',istep,x1,x2,z1,z2,b11,b12,b21,b22,bout))
      end
      fout = f12;
    end
    function fout = interpfield2(x,y,f,xq,yq)
      % INTERPFIELD3 Linear interpolation of 3D field.
      %   bout = interpfield(x,z,t,f,xq,yq,zq)
      %   x, z, t - grid
      %   f - field to interploate (defined on grid x,z,t)  
      %   xq, zq, tq - point to interpolate to, can be scalar values or vectors
      %
      nq = numel(xq);
      for iq = 1:nq
        % find two closest grid points in x and z
        %   11 ------- 21
        %    |  o      |      o is point (xq,zq)
        %    |         |
        %    |         |
        %   12---------22  

        % First do two 2D interpolations, then one 1D interpolation
        ix1 = find(x<xq(iq),1,'last');
        ix2 = find(x>xq(iq),1,'first');
        iy1 = find(y<yq(iq),1,'last');
        iy2 = find(y>yq(iq),1,'first');
        % pick out x and z values for these points
        x1 = x(ix1); % x
        x2 = x(ix2);
        y1 = y(iy1); % z
        y2 = y(iy2); 
        % pick out field values at the eight corners ()
        f111 = f(ix1,iy1,iz1);
        f112 = f(ix1,iy1,iz2);
        f121 = f(ix1,iy2,iz1);
        f122 = f(ix1,iy2,iz2);
        % distances between point 'o' and edges of box (x1,x2,z1,z2)
        xx = (xq(iq) - x1)/(x2-x1); %
        yy = (yq(iq) - y1)/(y2-y1);  
        zz = (zq(iq) - z1)/(z2-z1);
        % bilinear interpolation
        % you can solve for coefficients by matrix operations, but for now I just
        % used the expression directly: https://en.wikipedia.org/wiki/Bilinear_interpolation
        % Here are instructions for 3D: https://en.wikipedia.org/wiki/Trilinear_interpolation    
        % First do two 2D interpolations in yz-plane, then one 1D interpolation
        % along x (which order you do the planes in doesn't matter)
        % Bilinear
        f1 = f111*(1-xx)*(1-yy) + f121*xx*(1-yy) + f112*(1-xx)*yy + f122*xx*yy;

        % debug
        %disp(sprintf('istep = %g, [x1,x1,z1,z1] = [%g,%g,%g,%g],  [b11,b12,b21,b22] = [%g,%g,%g,%g], binterp = %g',istep,x1,x2,z1,z2,b11,b12,b21,b22,bout))
      end
      fout = f1;
    end
    function fout = interpfield1(x,f,xq)
      % INTERPFIELD3 Linear interpolation of 3D field.
      %   bout = interpfield(x,z,t,f,xq,yq,zq)
      %   x, z, t - grid
      %   f - field to interploate (defined on grid x,z,t)  
      %   xq, zq, tq - point to interpolate to, can be scalar values or vectors
      %
      nq = numel(xq);
      for iq = 1:nq
        % find two closest grid points in x and z
        %   11 ------- 21
        %    |  o      |      o is point (xq,zq)
        %    |         |
        %    |         |
        %   12---------22  

        % First do two 2D interpolations, then one 1D interpolation
        ix1 = find(x<xq(iq),1,'last');
        ix2 = find(x>xq(iq),1,'first');
        iy1 = find(y<yq(iq),1,'last');
        iy2 = find(y>yq(iq),1,'first');
        % pick out x and z values for these points
        x1 = x(ix1); % x
        x2 = x(ix2);
        y1 = y(iy1); % z
        y2 = y(iy2); 
        % pick out field values at the eight corners ()
        f111 = f(ix1,iy1,iz1);
        f112 = f(ix1,iy1,iz2);
        f121 = f(ix1,iy2,iz1);
        f122 = f(ix1,iy2,iz2);
        % distances between point 'o' and edges of box (x1,x2,z1,z2)
        xx = (xq(iq) - x1)/(x2-x1); %
        yy = (yq(iq) - y1)/(y2-y1);  
        zz = (zq(iq) - z1)/(z2-z1);
        % bilinear interpolation
        % you can solve for coefficients by matrix operations, but for now I just
        % used the expression directly: https://en.wikipedia.org/wiki/Bilinear_interpolation
        % Here are instructions for 3D: https://en.wikipedia.org/wiki/Trilinear_interpolation    
        % First do two 2D interpolations in yz-plane, then one 1D interpolation
        % along x (which order you do the planes in doesn't matter)
        % Bilinear
        f1 = f111*(1-xx)*(1-yy) + f121*xx*(1-yy) + f112*(1-xx)*yy + f122*xx*yy;

        % debug
        %disp(sprintf('istep = %g, [x1,x1,z1,z1] = [%g,%g,%g,%g],  [b11,b12,b21,b22] = [%g,%g,%g,%g], binterp = %g',istep,x1,x2,z1,z2,b11,b12,b21,b22,bout))
      end
      fout = f1;
    end  
  end
  methods (Access = protected)
    function out = get_binned_quantity(obj,field,iSpecies)
      % Get binned quantities from Smilei simulation
      % h5 structure: '/timestep00000000' (strange)
      
      % Check if field exists
      if not(any(contains(obj.fields,field)))
        error(sprintf('Unknown field ''%s''.',field))
      end
      % Find which file it is, should find one for each species
      iFile = find(contains(obj.attributes.deposited_quantity,field));      
      species = obj.attributes.deposited_species(iFile);
      iSp = obj.species{iSpecies};
      iFile = iFile(find(contains(species,iSp))); % file name starts at zero, so subtract one
      allFiles = dir([obj.particlebinning]);
      file = [allFiles(iFile).folder filesep allFiles(iFile).name];            
      
      infoPB = h5info(file);
      iterations = obj.iteration; % what is in the object (original, or might have been reduced by twci(pe)lim)
      nIter = numel(iterations);
      data = nan([obj.get_gridsize,nIter]);
      for iIter = 1:nIter
        iter = iterations(iIter);
        str_iter = sprintf('%08.0f',iter);
        data_tmp = h5read(file,...
          ['/timestep' str_iter],...
          [obj.grid{2}(1) obj.grid{1}(1)]',... % start indices
          [numel(obj.grid{2}) numel(obj.grid{1})]'); % number of counts
        data_tmp = data_tmp';
        data(:,:,iIter) = data_tmp;
      end
      out = data;
    end
    function out = get_field(obj,field)

% I want to make something that keeps data through loading sessions, for 
% example B that is used for many things.
%       persistent data_recycle
%       %data_recycle = [];
%       isSameObj = 0; isSameX = 0; isSameZ = 0; isSameT = 0;
%       data_to_recycle = {'bx','by','bz'};
%       data_to_recycle = {'','',''};
%       isFieldToRecycle = any(contains(data_to_recycle,field));
%       isNotEmptyRecycle = not(isempty(data_recycle));
%       if isNotEmptyRecycle && isfield(data_recycle,field)
%         isSameObj = strcmp(obj.file,data_recycle.(field).obj.file);
%         isSameX = isequal(obj.ix,data_recycle.(field).obj.ix);
%         isSameZ = isequal(obj.iz,data_recycle.(field).obj.iz);
%         isSameT = isequal(obj.it,data_recycle.(field).obj.it);
%       end
%       if all([isFieldToRecycle isNotEmptyRecycle isSameObj isSameX isSameZ isSameT isfield(data_recycle,field)])          
%         data = data_recycle.(field).data;        
%       else

        % get iterations
        iterations = obj.iteration;
        nIter = obj.length;
        % initialize matrix
        gridsize = obj.get_gridsize;
        data = nan([gridsize([3 2 1]),nIter]);
        data = nan([gridsize,nIter]);
        %disp(field)
        for iIter = 1:nIter
          %tt = tic;
          iter = iterations(iIter);
          str_iter = sprintf('%010.0f',iter);
          if strcmp(obj.software,'micPIC')
            data_tmp = h5read(obj.file,...
               ['/data/' str_iter '/' field],...
               [obj.grid{1}(1) obj.grid{2}(1)],... % start indices
               [numel(obj.grid{1}) numel(obj.grid{2})]); % number of counts
            %disp(sprintf('Reading %s: [%g %g] datapoints starting at [%g %g]',field,numel(obj.grid{1}),numel(obj.grid{2}),obj.grid{1}(1),obj.grid{2}(1)))
          elseif strcmp(obj.software,'Smilei') % Smilei
            data_tmp = h5read(obj.file,...
              ['/data/' str_iter '/' field],...
              [obj.grid{3}(1) obj.grid{2}(1) obj.grid{1}(1)]',... % start indices
              [numel(obj.grid{3}) numel(obj.grid{2}) numel(obj.grid{1})]'); % number of counts
            % dimension of data is saved as [nz ny nx], so permute it back here
            data_tmp = permute(data_tmp,[3 2 1 4]); % dimension of data is saved as 
            %  obj.info.Groups(1).Groups.obj.info.Groups(1).Groups.Datasets(1).DataspaceDatasets(1).Dataspace
          end
          data(:,:,:,iIter) = data_tmp;
          %toc(tt)
        end
%       end
      out = data;
%       if contains(field,data_to_recycle)
%         data_recycle.(field).obj = obj;
%         data_recycle.(field).data = out;
%       end      
    end
    function out = get_exp(obj,str_exp)
      % Loads data in expression and returns evaluated expression 
      % Example:
      % out  = obj.get_exp('log10(abs(Bx.*Bz))');
      eval_str = str_exp;
      list_methods = methods(obj);
      %list_fields = list_methods(cellfun(@(s) contains(str_exp,s),list_methods));
      list_fields = regexpi(str_exp,{'[a-z_A-Z]*'},'match');
      list_fields = intersect(list_fields{1},list_methods);

      varcount = 0;
      for ifield = 1:numel(list_fields)
        nchar = numel(list_fields{ifield});
        str_loc = strfind(str_exp,list_fields{ifield});
        for istr_ = 1:numel(str_loc)
          istr = str_loc(istr_);
          varcount = varcount + 1; 
          
          if numel(str_exp)>istr+nchar && not(isempty(strfind(str_exp(istr+nchar),'('))) && strfind(str_exp(istr+nchar),'(')
            ind1 = strfind(str_exp(istr+nchar:end),'(');
            ind2 = strfind(str_exp(istr+nchar:end),')');           
            indstr = str_exp(istr + nchar-1 + (ind1(1)+1:ind2(1)-1));            
            tmp_str = [list_fields{ifield} '(' indstr ')'];
            var = obj.(list_fields{ifield})(eval(indstr));
          else
            var = obj.(list_fields{ifield});
            tmp_str = list_fields{ifield};
          end 
          eval_str = strrep(eval_str,tmp_str,['allvars{' num2str(varcount) '}']);
          allvars{varcount} = var;
        end
      end
      out = eval(eval_str);           
    end  
  end
end