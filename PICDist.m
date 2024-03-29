classdef PICDist
  % Load PIC simulation data
  %   Does not contain all the data, but loads it in an easily accesible manner  
  %
  %   pic = PIC(h5FilePath)
  %   Bx = pic.Bx; % Bx is a (nt x nx x ny) matrix
  %   B = pic.B; % structure with 3 (nt x nx x ny) matrices  
  
  properties    
    nspecies
  end
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
    xe1_
    xe2_
    ze1_
    ze2_
    xi1_
    xi2_
    zi1_
    zi2_
    dxi_
    dzi_
    grid_
    indices_
    it_
    id_
    dists_
    tags_
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
    xe1
    xe2
    ze1
    ze2
    xi1
    xi2
    zi1
    zi2
    dxi
    dzi
    grid
    indices
    it
    id
    dists
    tags
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
  end
  
  methods
    function obj = PICDist(h5filePath)
      % sm = SMILEI(pathFields)
      % sm = SMILEI(pathFields,[],[]) - current implementation
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
      
      obj.file = h5filePath; 
      obj.info = h5info(h5filePath);
      obj.nspecies = h5readatt(obj.file, ['/'],'nSpecies');
      %obj.charge = obj.get_charge;
      %obj.mass = obj.get_mass;
      %uniqueMass = sort(unique(obj.mass));
      %obj.mime = uniqueMass(2)/uniqueMass(1); % second lightest/lightest
      %obj.teti = h5read(h5filePath,'/simulation_information/teti');
      %obj.wpewce = h5read(h5filePath,'/simulation_information/wpewce');
      
      obj.iteration = get_iterations(obj);
      obj.twpe = get_twpe(obj);
      obj.twci = get_twci(obj);
      obj.it = 1:numel(obj.iteration);
      obj.dists_ = obj.get_distlist;
      obj.tags_ = obj.get_tags;
      [xi1,xi2,zi1,zi2] = obj.get_locs;
      obj.xi1 = xi1;
      obj.xi2 = xi2;
      obj.zi1 = zi1;
      obj.zi2 = zi2;            
      
      for it = 1:numel(obj.it)
        obj.dxi{1,it} = obj.xi2{it} - obj.xi1{it};
        obj.dzi{1,it} = obj.zi2{it} - obj.zi1{it};
        %obj.dzi{iIter} = cell2mat(cellfun(@(x) diff(x),obj.zi{iIter},'UniformOutput',false));
      end
      
      
      obj.indices = cellfun(@(x) 1:numel(x),obj.dists,'UniformOutput',false);
      %obj.twpe = get_twpe(obj);      
      %obj.twci = obj.twpe/(obj.wpewce*obj.mime);
      %obj.indices_ = 1:numel(obj.iteration);
      
      %obj.fields_ = get_fields(obj);
      %obj.xe = h5read(h5filePath,'/simulation_information/xe'); % de
      %obj.ze = h5read(h5filePath,'/simulation_information/ze'); % de
      %obj.xi = obj.xe/sqrt(obj.mime);
      %obj.zi = obj.ze/sqrt(obj.mime);
      %obj.grid = {1:1:numel(obj.xe),1:1:numel(obj.ze)}; % originally, complete grid
      %obj.ix = 1:1:numel(obj.xe);
      %obj.iz = 1:1:numel(obj.ze);
      %obj.it = 1:1:numel(obj.twpe);
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
          obj.iteration_ = builtin('subsref',obj.iteration_,s);
          %obj.twpe_ = builtin('subsref',obj.twpe,s);
          %obj.twci_ = builtin('subsref',obj.twci,s);
          obj.indices_ = builtin('subsref',obj.indices_,s);
          obj.dists_ = builtin('subsref',obj.dists_,s);
          obj.xi1_ = builtin('subsref',obj.xi1_,s);
          obj.xi2_ = builtin('subsref',obj.xi2_,s);
          obj.zi1_ = builtin('subsref',obj.zi1_,s);
          obj.zi2_ = builtin('subsref',obj.zi2_,s);
          obj.dxi_ = builtin('subsref',obj.dxi_,s);
          obj.dzi_ = builtin('subsref',obj.dzi_,s);
          obj.it_ = builtin('subsref',obj.it_,s);
          obj.twpe_ = builtin('subsref',obj.twpe_,s);
          obj.twci_ = builtin('subsref',obj.twci_,s);
         
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
    function obj = find(obj,field,expr,value)
      % NOT IMPLEMENTED
      for it = obj.it
        [I,J] = find(eval(['obj.' field expr num2str(value)]))
        
        
      end      
    end
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
        x1 = obj.xi1{it};
        x2 = obj.xi2{it};
        if strcmp(method,'center')
          xtmp = (x1+x2)/2;
        end
        i1 = find(xtmp >= value(1)); 
        i2 = find(xtmp <= value(2)); 
        inds{it} = intersect(i1,i2);
      end            
      
      % Update grid and indices
      obj = obj.update_inds(inds);
    end
    function obj = zlim(obj,value,varargin)
      % Get subset of z
      nt = obj.nt;
      nd = obj.nd;
      method = 'center';
      
      for it = 1:nt
        z1 = obj.zi1{it};
        z2 = obj.zi2{it};
        if strcmp(method,'center')
          ztmp = (z1+z2)/2;
        end
        i1 = find(ztmp >= value(1)); 
        i2 = find(ztmp <= value(2)); 
        inds{it} = intersect(i1,i2);
      end            
      
      % Update grid and indices
      obj = obj.update_inds(inds);
    end
    function obj = xfind(obj,value)
      % pick distributions with centers corresponding to exact value           
      
      for it = 1:obj.nt
        center = (obj.xi1{it} + obj.xi2{it})/2;
%         for ii = 1:numel(center)
%           if find(value==center(ii),1,'first')
%             IA(ii) = 1;
%           else
%             IA(ii) = 0;
%           end
%         end
        [IA,IB] = ismembertol(center,value,1e-5);
        inds{it} = find(IA);
      end
      obj = obj.update_inds(inds);
    end
    function obj = zfind(obj,value)
      % pick distributions with centers corresponding to exact value           
      
      for it = 1:obj.nt        
        center = (obj.zi1{it} + obj.zi2{it})/2;
        [IA,IB] = ismembertol(center,value,1e-3);
        inds{it} = find(IA);
      end      
      obj = obj.update_inds(inds);
    end
    function obj = findx(obj,value)
      obj = xfind(obj,value);
    end
    function obj = findz(obj,value)
      obj = zfind(obj,value);
    end
    function obj = findtag(obj,value)
      % Finds and picks special tags
      tags = obj.tags;
      iVals = numel(value);
      for it = 1:obj.nt
        isTag = strcmp(tags{it},value);
        inds{it} = find(isTag);
      end      
      obj = obj.update_inds(inds);
    end
    function obj = dxlim(obj,value)
      % Get subset of dx (box size)
      nt = obj.nt;
      nd = obj.nd;      
      
      for it = 1:nt
        dtmp = obj.dxi{it};        
        i1 = find(dtmp >= value(1)); 
        i2 = find(dtmp <= value(2)); 
        inds{it} = intersect(i1,i2);
      end            
      
      % Update grid and indices
      obj = obj.update_inds(inds);   
    end
    function obj = dzlim(obj,value)
      % Get subset of dz (box size)
      nt = obj.nt;
      nd = obj.nd;      
      
      for it = 1:nt
        dtmp = obj.dzi{it};        
        i1 = find(dtmp >= value(1)); 
        i2 = find(dtmp <= value(2)); 
        inds{it} = intersect(i1,i2);
      end            
      
      % Update grid and indices
      obj = obj.update_inds(inds);   
      
    end
    function obj = twpelim(obj,value,varargin)
      % Get subset of twpe
      inds = obj.ind_from_lim(obj.twpe_,value,varargin{:});
      obj = obj.update_inds_time(inds);
    end
    function obj = twcilim(obj,value,varargin)
      % Get subset of twci
      inds = obj.ind_from_lim(obj.twci,value,varargin{:});
      obj = obj.update_inds_time(inds);
    end
    function obj = twcifind(obj,value)
      % pick distributions with centers corresponding to exact value           
      
      inds = find(ismember(obj.twci,value));
      obj = obj.update_inds_time(inds);
    end
    
    function obj = update_inds(obj,inds)
      % xi
      % zi
      % dists
      % indices
      for it = 1:obj.nt
        obj.xi1_{it} = obj.xi1_{it}(inds{it});
        obj.xi2_{it} = obj.xi2_{it}(inds{it});
        obj.zi1_{it} = obj.zi1_{it}(inds{it});
        obj.zi2_{it} = obj.zi2_{it}(inds{it});
        obj.dxi_{it} = obj.dxi_{it}(inds{it});
        obj.dzi_{it} = obj.dzi_{it}(inds{it});
                
        obj.dists_{it} = obj.dists_{it}(inds{it});
        obj.indices_{it} = obj.indices_{it}(inds{it});
        obj.tags_{it} = obj.tags_{it}(inds{it});        
      end
      % run through and update time inds if any time is empty: empty times
      % are removed
      emptytimes = cellfun(@(x)isempty(x),obj.xi1,'UniformOutput',false);
      time_inds = find(not([emptytimes{:}]));
      obj = obj.update_inds_time(time_inds);      
    end
    function obj = update_inds_time(obj,inds)
      obj.twpe_ = obj.twpe_(inds);
      obj.twci_ = obj.twci_(inds);
      obj.it_ = obj.it_(inds);
      obj.indices_ = obj.indices_(inds);
      obj.iteration_ = obj.iteration_(inds);
      obj.xi1 = obj.xi1(inds);
      obj.xi2 = obj.xi2(inds);
      obj.zi1 = obj.zi1(inds);
      obj.zi2 = obj.zi2(inds);
      obj.dxi = obj.dxi(inds);
      obj.dzi = obj.dzi(inds);
      obj.dists = obj.dists(inds); 
      obj.tags = obj.tags(inds);      
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
    function varargout = plot_map(obj,varargin)
      % PICDist.plot_map Plots a map of all distributions in new figure.
      %   h = plot_map(obj,iSpecies,sumdim,options);
      %
      %   Options:
      %     'bline',PICobj - draw in plane magnetic field line
      %     'v',PICobj - plot in plane velocity moment
      %     'exb',PICobj - plot in plane ExB velocity
      %     'log' - plot log10(data)
      %     'diff',ispecies - plot difference between phase space density
      %           of two species, e.g. obj.plot_map(ispecies,sumdim,diffSpecies) plots f_ispecies-f_diffspecies (or maybe f_diffspecies-f_ispecies, verify)
      %     'ratio',ispecies - same as 'diff', but just divide them two
      %           instead, e.g. obj.plot_map(ispecies,sumdim,diffSpecies) 
      %           plots f_ispecies/f_diffspecies, so obj.plot_map(1,sumdim,[1 3])
      %           is f1/(f1+f3)
      %
      %   Output:
      %     h.ax - handle to all axes
      %     h.leg - handle to all text labels (showing location)
      %     h.links - x and y axes are linked (if you zoom on one you zoom 
      %               on all)
      %   
      %
      
      doLine = 0;
      doLineDrift = 0;
      doDot = 1;
      strDot = {};
      doColorbar = 1;
      doV = 0;
      fontsize = 7;
      border = 0;
      %doBDir = 1;
      doLog = 0;
      ticks = -15:1:15;
      doNumber = 1;
      doDiff = 0;
      doFrac = 0;
      doRatio = 0;
      doForce = 0;
      doCurvV = 0;
      cmap = pic_colors('candy4');
      doSmooth = 0;
      plotInAxes = 0;
      doLabel = 1;
      doNaN = 0;
      doExB = 0;
      doCS_fieldaligned = 0;
      doVectors = 0;
      doPIntegrand = 0;
      doContour = 0;
      
      [ax,args,nargs] = irf.axescheck(varargin{:});       
      iSpecies = args{1}; args = args(2:end); nargs = nargs - 1;
      sumdim = args{1}; args = args(2:end); nargs = nargs - 1;
                  
      if not(isempty(ax)); plotInAxes = 1; doCompact = 0; end
      
      % check if input PICDist obj has a single time, otherwise take first
      % time, and warn
      if obj.nt > 1
        warning('PICDist has more than one times, selecting the first: twpe = %g, twci = %g',obj.twpe(1),obj.twci(1))
        obj = obj(1);
      end
      
      have_options = 0;
      if not(isempty(args))
        have_options = 1;
      end
      while have_options
        l = 1;
        switch lower(args{1})
          case 'bline' % plot a line defined by argument
            doLine = 1;
            l = 2;
            pic = args{2};            
          case 'line' % plot a line defined by argument
            doLine = 1;
            l = 2;
            line = args{2};
          case 'line_drift' % shift the direction line by some velocity value
            doLine = 1;
            doLineDrift = 1;
            l = 3;
            line = args{2};
            line_shift = args{3};
          case 'v'
            doV = 1;
            pic = args{2};
            l = 2;
          case 'exb'
            doExB = 1;
            pic = args{2};
            l = 2;            
          case 'log'
            doLog = 1;            
            l = 1;
          case 'diff' %
            doDiff = 1;
            iSpeciesDiff = args{2};
            l = 2;
          case 'frac' %
            doFrac = 1;
            iSpeciesDiff = args{2};
            l = 2;
          case 'ratio' %
            doRatio = 1;
            iSpeciesDiff = args{2};
            l = 2;
          case 'force'
            doForce = 1;
            force_exp = args{2};
            pic = args{3};            
            l = 3;
            cmap = pic_colors('blue_red');
          case 'smooth'
            doSmooth = 1;
            npSmooth = args{2};
            l = 2;
          case 'curv'
            doCurvV = 1;
            pic = args{2}{1};
            kCurv = args{2}{2};
            l = 2;
          case 'nolabel'
            doLabel = 0;
            l = 1;
          case 'nan'
            doNaN = 1;
            l = 1;
          case 'fieldaligned_cs'
            doCS_fieldaligned = 1;
            pic = args{2};
            l = 2;            
          case 'vectors'
            doVectors = 1;
            vectors = args{2};
            pic = args{3};
            l = 3;
          case 'off-diag'
            doPIntegrand = 1;
            l = 2;
            mass = args{2};
          case 'p-integrand'
            doPIntegrand = 1;
            l = 2;  
            P_comps = args{2};
            cmap = pic_colors('blue_red');
        end
        args = args((1+l):end);
        if isempty(args); break; end
      end
      
      if mod(iSpecies,2) == 1 % ions are uneven number
        ticks = [-20:1:20];
      else
        ticks = [-20:1:20];
      end
      strspecies = '[';
      for isp = 1:numel(iSpecies)
        strspecies = [strspecies num2str(iSpecies(isp)) ','];
      end
      strspecies(end) =']';
            
      if not(plotInAxes)      
        xlim = [min([obj.xi1_{1}(:)]) max([obj.xi2_{1}(:)])];
        zlim = [min([obj.zi1_{1}(:)]) max([obj.zi2_{1}(:)])];      
   %     fig_position = get(0,'screensize'); %[1 330 2555 1015]; 
   %     fig = figure;
   %     fig.Position = fig_position; 
        fig = gcf;
        border = [0.05,0.1,0.07,0.05]; % left, bottom, right, top
      end
      idist = 0;
      %tic
      h = gobjects(0);
      hleg = gobjects(0);
      for id = obj.indices{1}                        
        idist = idist + 1;      

        % Load data                  
        f = obj.fxyz(1,idist,iSpecies,sumdim); % sum over 3rd dim
        if doDiff
          fdiff = obj.fxyz(1,idist,iSpeciesDiff,sumdim); % sum over 3rd dim
          %f.f = f.f./fdiff.f;
          f.f = f.f - fdiff.f;
          f.f(f.f==0) = NaN;
          % This becomes wrong because log10 of something positive can be
          % negative
          %if doLog
          %  f.f = sign(f.f).*log10(abs(f.f));
          %end
          cmap = pic_colors('blue_red');
        end
        if doRatio
          fdiff = obj.fxyz(1,idist,iSpeciesDiff,sumdim); % sum over 3rd dim
          %f.f = f.f./fdiff.f;
          f.f = f.f./fdiff.f;
          f.f(fdiff.f == 0) = NaN;
          cmap = pic_colors('blue_red');
        end
        if doFrac
          fdiff = obj.fxyz(1,idist,iSpeciesDiff,sumdim); % sum over 3rd dim
          %f.f = f.f./fdiff.f;
          if doLog
            f.f = log10(f.f) - log10(fdiff.f);
          else
            f.f = f.f./fdiff.f;
          end
          cmap = pic_colors('blue_red');
        end
        if doForce
          switch force_exp
             case 'EvBy'
               Ey = pic.get_exp('Ey');
            case 'vx*Bz'
              pic_tmp = pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z);
              Bz = mean(mean(pic_tmp.Bz));
              [V1,V2] = ndgrid(f.v,f.v);
              % assume V1 is vx
              force = V1*Bz;
            case 'vx*Bz-Ey'
              pic_tmp = pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z);
              Bz = mean(mean(pic_tmp.Bz));
              Ey = mean(mean(pic_tmp.Ey));
              [V1,V2] = ndgrid(f.v,f.v);
              % assume V1 is vx
              force = V1*Bz-Ey;
            case '-vy*Bz-Ex'
              pic_tmp = pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z);
              Bz = mean(mean(pic_tmp.Bz));
              Ex = mean(mean(pic_tmp.Ex));
              [V1,V2] = ndgrid(f.v,f.v);
              % assume V1 is vx
              force = -V2*Bz-Ex;
          end

        end
        if doPIntegrand
          [V1,V2] = ndgrid(f.v,f.v);  
          n = sum(sum(f.f))*f.dv^2;
          fv1 = f.f.*V1;
          fv2 = f.f.*V2;
          vbulk1 = sum(fv1(:))*(f.dv^2)/n;
          vbulk2 = sum(fv2(:))*(f.dv^2)/n;
          f.f = f.f.*(V1-vbulk1).*(V2-vbulk2);
          colormap(pic_colors('blue_red'))
        end
        %if not(all(xlo>xlim(1) && xhi<xlim(2) && zlo>zlim(1) && zhi<zlim(2)))
        %  disp([sprintf('%.2f %.2f %.2f %.2f outside of box',xlo,xhi,zlo,zhi)])
        %  continue
        %end
%         xloc = (f.x(1)+f.x(2))/2;
%         zloc = (f.z(1)+f.z(2))/2;
%         dx = f.x(2)-f.x(1);
%         dz = f.z(2)-f.z(1);
        %disp(sprintf('%.2f ',f.x(1),f.x(2),f.z(1),f.z(2)))
        if not(plotInAxes)
          axes_position = [(f.x(1)-xlim(1))/diff(xlim) ...
                           (f.z(1)-zlim(1))/diff(zlim) ...
                           (f.x(2)-f.x(1))/diff(xlim) ...
                           (f.z(2)-f.z(1))/diff(zlim)];
          axes_position(1) = axes_position(1)*(1-border(1)-border(3)) + border(1);
          axes_position(2) = axes_position(2)*(1-border(2)-border(4)) + border(2);
          axes_position(3) = axes_position(3)*(1-border(1)-border(3));
          axes_position(4) = axes_position(4)*(1-border(2)-border(4));
        end
        
        pause(0.1)  
        switch sumdim
          case 1
            xaxstr = 'y';
            yaxstr = 'z';
          case 2
            xaxstr = 'x';
            yaxstr = 'z';
          case 3
            xaxstr = 'x';
            yaxstr = 'y';
        end

        if 1 % plane
          nrows = 1;
          ncols = 3;
          npanels = nrows*ncols;
          
          if not(plotInAxes)
            figure(fig);          
            hca = axes('Position',axes_position);
            h(end+1) = hca;
            hca = gca;
          else
            hca = ax(idist);
            h(idist) = hca;
          end
          %if doForce
%             switch force_exp
%               case 'EvBx', if not(sumdim == 1), continue, end
%                 pic_tmp = pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z);
%                 Ex = squeeze(mean(mean(pic_tmp.get_exp('Ey'))));
%                 By = squeeze(mean(mean(pic_tmp.get_exp('Bx'))));
%                 Bz = squeeze(mean(mean(pic_tmp.get_exp('Bz'))));
%                 [VY,VZ] = meshgrid(f.v,f.v);
%                 force = Ex + VZ*By - VY*Bz;
%               case 'EvBy', if not(sumdim == 2), continue, end
%                 pic_tmp = pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z);
%                 Ey = squeeze(mean(mean(pic_tmp.get_exp('Ey'))));
%                 Bx = squeeze(mean(mean(pic_tmp.get_exp('Bx'))));
%                 Bz = squeeze(mean(mean(pic_tmp.get_exp('Bz'))));
%                 [VX,VZ] = meshgrid(f.v,f.v);
%                 force = Ey + VZ*Bx - VX*Bz;
%               case 'EvBz', if not(sumdim == 3), continue, end
%                 pic_tmp = pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z);
%                 Ez = squeeze(mean(mean(pic_tmp.get_exp('Ey'))));
%                 Bx = squeeze(mean(mean(pic_tmp.get_exp('Bx'))));
%                 By = squeeze(mean(mean(pic_tmp.get_exp('By'))));
%                 [VX,VY] = meshgrid(f.v,f.v);
%                 force = Ez + VX*By - VY*Bx;
%               case 'vBy', if not(sumdim == 2), continue, end
%                 pic_tmp = pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z);               
%                 Bx = squeeze(mean(mean(pic_tmp.get_exp('Bx'))));
%                 Bz = squeeze(mean(mean(pic_tmp.get_exp('Bz'))));
%                 [VX,VZ] = meshgrid(f.v,f.v);
%                 force = VZ*Bx - VX*Bz;
%               case 'Ey'
%                 pic_tmp = pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z);
%                 Ey = squeeze(mean(mean(pic_tmp.get_exp('Ey'))));
%                 [VX,VZ] = meshgrid(f.v,f.v);
%                 force = Ey;
%               otherwise
%                 [VX,VZ] = meshgrid(f.v,f.v);
%                 force = VX*0;
%             end
%             force(f.f==0) = nan;
%             imagesc(hca,f.v,f.v,force')
          if doLog
            title_str = sprintf('log_{10}f(v_%s,v_%s,twci=%g), species = %s',xaxstr,yaxstr,obj.twci,strspecies);
          else
            title_str = sprintf('f(v_%s,v_%s), species = %s',xaxstr,yaxstr,strspecies);
          end
            
          if doLog && any([doFrac doRatio])
            fplot = log10(f.f);            
            cmap = pic_colors('blue_red');
            %vsr_str = 
          elseif doLog && doDiff    
            fplot = sign(f.f).*log10(abs(f.f));
          elseif doLog
            fplot = log10(f.f);
          else
            fplot = f.f;
          end          
          if doSmooth
            fplot = smooth2(fplot,npSmooth);
          end
          if doNaN
            fplot(fplot==0) = NaN;
          end
          if doForce
            fplot = force;
            doContour = 1;
          end
          imagesc(hca,f.v,f.v,fplot')
          
          if doContour
            hold(hca,'on')
            contour(hca,f.v,f.v,smooth2(f.f,3)','k')
            hold(hca,'off')
          end
          
          %hca.XLabel.String = '';
          %hca.YLabel.String = '';
          hca.YDir = 'normal';        
          %hca.XLim = vlim(ispecies)*[-1 1];
          %hca.YLim = vlim(ispecies)*[-1 1];
          hca.XGrid = 'on';
          hca.YGrid = 'on';
          hca.XTick = ticks;
          hca.YTick = ticks;              
          %hca.XTickLabel = [];
          %hca.YTickLabel = [];
          hca.Box = 'on';
          colormap(hca,cmap)
          %irf_legend(hca,{sprintf('x=%.1f, z=%.1f',xloc,zloc);sprintf('B=[%.2f,%.2f,%.2f]',Bloc.x,Bloc.y,Bloc.z)},[0.01 0.99],'color',[0 0 0],'fontsize',9)
          if doLabel
            hleg_ = irf_legend(hca,{...
              sprintf('x=[%.2f,%.2f]',f.x(1),f.x(2));...
              sprintf('z=[%.2f,%.2f]',f.z(1),f.z(2))},...
              [0.01 0.01],'color',[0 0 0],'fontsize',fontsize);
            hleg(size(hleg,1)+1,:) = hleg_;
          end
          if doPIntegrand % make caxes symmetric
            hca.CLim = max(abs(hca.CLim))*[-1 1];
            sumf = sum(fplot(:));
            if sumf >= 0 
              sumf_color = 'r';
            else
              sumf_color = 'b';
            end
            irf_legend(hca,{sprintf('sum = %g',sumf*mass)},[0.02 0.98],'color',sumf_color,'fontsize',fontsize)
          end
          if not(plotInAxes) % xy-labels
            if axes_position(1) == border(1)
              hca.YLabel.String = sprintf('v_{%s}',yaxstr);            
            else
              hca.YLabel.String = '';
              hca.YTickLabel = [];
            end
            if axes_position(2) == border(2)
              hca.XLabel.String = sprintf('v_{%s}',xaxstr);
            else
              hca.XLabel.String = '';
              hca.XTickLabel = [];
            end
            if axes_position(1) == border(1) && abs(axes_position(2)+axes_position(4) - (1-border(4))) < 1e-5            
            if doColorbar
              doColorbar = 0; 
              hb = colorbar('peer',hca);
              barwidth = 0.01;
              hb.Position = [1-border(3)+0.5*barwidth border(2) barwidth 1-border(2)-border(4)];
              hb.YLabel.String = title_str;
            else
              hca.Title.String = title_str;
            end                        
            end
          else
            hca.YLabel.String = sprintf('v_{%s}',yaxstr);
            hca.XLabel.String = sprintf('v_{%s}',xaxstr);
          end
          if doLine
            pic_tmp = pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z);
            Bx = mean(mean(pic_tmp.Bx));
            By = mean(mean(pic_tmp.By));
            Bz = mean(mean(pic_tmp.Bz));
            Babs = sqrt(Bx.^2+By.^2+Bz.^2);
            bx = Bx/Babs;
            by = By/Babs;
            bz = Bz/Babs;
            switch sumdim
              case 2 % xz
                line_slope = (bz/bx);
              case 1 % yz
                line_slope = (bz/by);
              case 3 % xy
                line_slope = (by/bx);  
            end                                    
            xx = min(hca.XLim(2)*[1 1/abs(line_slope)])*[-1 1];
            hold(hca,'on')                
            hBline = plot(hca,xx,xx*line_slope,'linewidth',0.5,'color',[0.5 0.5 0.5]);
            hold(hca,'off')
          end
          if doV
            pic_tmp = pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z);
            switch sumdim
              case 2 % xz
                v1 = mean(mean(pic_tmp.vx(iSpecies)));                
                v2 = mean(mean(pic_tmp.vz(iSpecies)));
              case 1 % yz
                v1 = mean(mean(pic_tmp.vy(iSpecies)));                
                v2 = mean(mean(pic_tmp.vz(iSpecies)));
              case 3 % xy
                v1 = mean(mean(pic_tmp.vx(iSpecies)));                
                v2 = mean(mean(pic_tmp.vy(iSpecies)));
            end
            hold(hca,'on')                
            %hVDot = plot(hca,v1,v2,'Marker','x','color',[0 0 0]);
            hVDot = scatter(hca,v1,v2,30,[0 0 0],'Marker','x');
            hold(hca,'off')
          end
          if doExB
            pic_tmp = pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z);
%             switch sumdim
%               case 2 % xz
%                 v1 = mean(mean(pic_tmp.vExBx));                
%                 v2 = mean(mean(pic_tmp.vExBz));
%               case 1 % yz
%                 v1 = mean(mean(pic_tmp.vExBy));                
%                 v2 = mean(mean(pic_tmp.vExBz));
%               case 3 % xy
%                 v1 = mean(mean(pic_tmp.vExBx));                
%                 v2 = mean(mean(pic_tmp.vExBy));
%             end
            vExB1 = mean(mean(pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z).(['vExB' xaxstr])));
            vExB2 = mean(mean(pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z).(['vExB' yaxstr])));            
            hold(hca,'on')                
            %hVDot = plot(hca,v1,v2,'Marker','x','color',[0 0 0]);
            hVExB = scatter(hca,vExB1,vExB2,20,[0 0 0],'Marker','o');
            hold(hca,'off')
          end
          if doNumber
            
          end
          if doForce
            clim = hca.CLim;
            hold(hca,'on')
            switch force_exp
              case 'EvBx', if not(sumdim == 1), continue, end
                pic_tmp = pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z);
                Ex = squeeze(mean(mean(pic_tmp.get_exp('Ey'))));
                By = squeeze(mean(mean(pic_tmp.get_exp('Bx'))));
                Bz = squeeze(mean(mean(pic_tmp.get_exp('Bz'))));
                [VY,VZ] = meshgrid(f.v,f.v);
                force = Ex + VZ*By - VY*Bz;
              case 'EvBy', if not(sumdim == 2), continue, end
                pic_tmp = pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z);
                Ey = squeeze(mean(mean(pic_tmp.get_exp('Ey'))));
                Bx = squeeze(mean(mean(pic_tmp.get_exp('Bx'))));
                Bz = squeeze(mean(mean(pic_tmp.get_exp('Bz'))));
                [VX,VZ] = meshgrid(f.v,f.v);
                force = Ey + VZ*Bx - VX*Bz;
              case 'EvBz', if not(sumdim == 3), continue, end
                pic_tmp = pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z);
                Ez = squeeze(mean(mean(pic_tmp.get_exp('Ey'))));
                Bx = squeeze(mean(mean(pic_tmp.get_exp('Bx'))));
                By = squeeze(mean(mean(pic_tmp.get_exp('By'))));
                [VX,VY] = meshgrid(f.v,f.v);
                force = Ez + VX*By - VY*Bx;
              case 'vBy', if not(sumdim == 2), continue, end
                pic_tmp = pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z);               
                Bx = squeeze(mean(mean(pic_tmp.get_exp('Bx'))));
                Bz = squeeze(mean(mean(pic_tmp.get_exp('Bz'))));
                [VX,VZ] = meshgrid(f.v,f.v);
                force = VZ*Bx - VX*Bz;
              case 'Ey'
                pic_tmp = pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z);
                Ey = squeeze(mean(mean(pic_tmp.get_exp('Ey'))));
                [VX,VZ] = meshgrid(f.v,f.v);
                force = Ey;
              otherwise
                [VX,VZ] = meshgrid(f.v,f.v);
                force = VX*0*nan;
            end
            %force(f.f==0) = nan;
            contour(hca,f.v,f.v,force',[0 0],'k')
            hold(hca,'off')
            hca.CLim = clim;
          end
          if doCurvV
            hold(hca,'on')
            Babs = mean(mean(pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z).Babs));
            curvbrad = mean(mean(pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z).curvbrad));
            vExB1 = mean(mean(pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z).(['vExB' xaxstr])));
            vExB2 = mean(mean(pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z).(['vExB' yaxstr])));            
            for iCurvV = 1:numel(kCurv)
              v1 = vExB1+curvbrad*Babs*kCurv(iCurvV)*cosd(0:360);
              v2 = vExB2+curvbrad*Babs*kCurv(iCurvV)*sind(0:360);
              hCurvV = plot(hca,v1,v2,'color',[0 0 0]);
              scatter(hca,vExB1,vExB2,40,[0 0 0],'marker','.');
              %plot(hca,vExB1,vExB2,'marker','*','color',[0 0 0],'markersize',3);
            end
            hold(hca,'off')            
          end
          
          if doCS_fieldaligned
            Bx = mean(mean(pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z).Bx));
            By = mean(mean(pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z).By));
            Bz = mean(mean(pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z).Bz));
            Ex = mean(mean(pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z).Ex));
            Ey = mean(mean(pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z).Ey));
            Ez = mean(mean(pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z).Ez));
            
            b = [Bx,By,Bz]/sqrt(Bx.^2+By.^2+Bz.^2);
            e = [Ex,Ey,Ez]/sqrt(Ex.^2+Ey.^2+Ez.^2);
            exb = cross(e,b); exb = exb/norm(exb);
            bxexb = cross(b,exb); bxexb = bxexb/norm(bxexb);
            
            s = hca.XLim(2)*0.7;            
                        
            
            hold(hca,'on')
            switch sumdim
              case 1 % f(vy,vz)
               quiver(hca,0,0,b(2)*s,b(3)*s,0,'k')
               quiver(hca,0,0,exb(2)*s,exb(3)*s,0,'k')
               quiver(hca,0,0,bxexb(2)*s,bxexb(3)*s,0,'k')   
               
               text(hca,b(2)*s,b(3)*s,sprintf('b = [%.2f,%.2f,%.2f]',b(1),b(2),b(3)))
               text(hca,exb(2)*s,exb(3)*s,sprintf('exb = [%.2f,%.2f,%.2f]',exb(1),exb(2),exb(3)))
               text(hca,bxexb(2)*s,bxexb(3)*s,sprintf('bx(exb) = [%.2f,%.2f,%.2f]',b(1),b(2),b(3)))
              case 2 % f(vx,vz)
               quiver(hca,0,0,b(1)*s,b(3)*s,0,'k')
               quiver(hca,0,0,exb(1)*s,exb(3)*s,0,'k')
               quiver(hca,0,0,bxexb(1)*s,bxexb(3)*s,0,'k')   
               
               text(hca,b(1)*s,b(3)*s,sprintf('B = [%.2f,%.2f,%.2f]',b(1),b(2),b(3)))
               text(hca,exb(1)*s,exb(3)*s,sprintf('exb = [%.2f,%.2f,%.2f]',b(1),b(2),b(3)))
               text(hca,bxexb(1)*s,bxexb(3)*s,sprintf('bx(exb) = [%.2f,%.2f,%.2f]',b(1),b(2),b(3)))
              case 3 % f(vx,vy)
               quiver(hca,0,0,b(1)*s,b(2)*s,0,'k')
               quiver(hca,0,0,exb(1)*s,exb(2)*s,0,'k')
               quiver(hca,0,0,bxexb(1)*s,bxexb(2)*s,0,'k')   
               
               text(hca,b(1)*s,b(2)*s,sprintf('B = [%.2f,%.2f,%.2f]',b(1),b(2),b(3)))
               text(hca,exb(1)*s,exb(2)*s,sprintf('exb = [%.2f,%.2f,%.2f]',b(1),b(2),b(3)))
               text(hca,bxexb(1)*s,bxexb(2)*s,sprintf('bx(exb) = [%.2f,%.2f,%.2f]',b(1),b(2),b(3)))
            end                        
            hold(hca,'off')
          end

          if doVectors
            pic_tmp = pic.twpelim(obj.twpe,'exact').xlim(f.x).zlim(f.z);
            for iVec = 1:numel(vectors)
              
            end
          end
          %print('-dpng','-r200',[savedir_root sub_dir '/' strprint '.png']);
          drawnow
          %pause(1)
        end
      end
      if doColorbar && not(plotInAxes)
        hb = colorbar('peer',hca);
        barwidth = 0.01;
        hb.Position = [1-border(3)+0.5*barwidth border(2) barwidth 1-border(2)-border(4)];
        hb.YLabel.String = title_str;
      end
        %toc      
      hlinks = linkprop(h,{'XLim','YLim','XTick','YTick'});
      varargout{1}.ax = h;
      varargout{1}.leg = hleg;
      varargout{1}.links = hlinks;
      
    end
    function varargout = plot_map_pitch(obj,iSpecies,pic,varargin)
      
      doSmooth = 0;
      
      have_options = 0;
      if not(isempty(varargin))
        args = varargin;
        have_options = 1;
      end
      while have_options
        l = 1;
        switch lower(args{1})
          case 'smooth'
            doSmooth = 1;
            npSmooth = args{2};
            l = 2;
        end
        args = args((1+l):end);        
        if isempty(args); break; end
      end
      
      
    end
    function varargout = plot_boxes(obj,varargin)
      % PDIST.PLOT_BOXES Plot boxes of distributions.
      
      doDe = 0;
      doTagColor = 0;
      color = [0 0 0];
      doPlotIntoAxes = 0;
      [ax,args,nargs] = axescheck(varargin{:});
            
      have_options = 0;
      if not(isempty(args))
        have_options = 1;
      end
      while have_options
        l = 1;
        switch lower(args{1})
          case 'de'
            doDe = 1;
            mime = args{2};
            l = 2;
          case 'color' % plot a line defined by argument
            color = args{2};
            if ischar(color) && strcmp(color,'tags')
              doTagColor = 1;
            end
            l = 2;
        end
        args = args((1+l):end);
        if isempty(args); break; end
      end
      % Get min and max values in order to set common XLim and YLim
      nt = obj.nt;
      
      if not(nt == numel(ax))
        warning('Numer of timesteps not equal to number of axes. Plotting into new figure.')
        figure;
      else
        doPlotIntoAxes = 1;
      end
      
      
      x1all = []; x2all = []; z1all = []; z2all = [];
      for it = 1:nt
        x1all = [x1all obj.xi1{it}];
        x2all = [x2all obj.xi2{it}];
        z1all = [z1all obj.zi1{it}];
        z2all = [z2all obj.zi2{it}];
      end
      xmin = min(x1all); 
      xmax = max(x2all);
      zmin = min(z1all); 
      zmax = max(z2all);
      
      for it = 1:nt
        ihsave = [];
        skipTagColor = 0;
        if doTagColor % set up colors based on tags
          tags = obj.tags{it};
          emptyCell = find(cellfun(@(s) isempty(s), tags));
          tags(emptyCell) = repmat({' '},numel(emptyCell),1);
          if isempty([tags{:}])
            color = [0 0 0];
            skipTagColor = 1;
          else
            uniquetags = unique(tags);
            ntags = numel(uniquetags);
            cmap = pic_colors('pasteljet');
            colors = interp1((0:size(cmap,1)-1)/size(cmap,1),cmap,(0:ntags-1)/ntags);          
          end
        end
        if doPlotIntoAxes
          hca = ax(it);
          hold(hca,'on')
        else
          hca = subplot(nt,1,it); ax(it) = hca;
        end
        for idist = 1:numel(obj.xi1{it})
          xplot = [obj.xi1{it}(idist) obj.xi2{it}(idist) obj.xi2{it}(idist) obj.xi1{it}(idist) obj.xi1{it}(idist)];
          zplot = [obj.zi1{it}(idist) obj.zi1{it}(idist) obj.zi2{it}(idist) obj.zi2{it}(idist) obj.zi1{it}(idist)];
          if doTagColor && not(skipTagColor) % set up colors
            tag = tags{idist};
            iuniquetag = find(cellfun(@(s) ~isempty(strfind(tag, s)), uniquetags));
            color = colors(iuniquetag,:);        
          end
          if doDe
            hbox(idist) = plot(hca,xplot*sqrt(mime),zplot*sqrt(mime),'color',color,'linewidth',1);
          else
            hbox(idist) = plot(hca,xplot,zplot,'color',color,'linewidth',1);
          end
          if doTagColor && not(skipTagColor)% Save first handle of each tag            
            ihsave(iuniquetag) = idist;
          end
          if idist == 1; hold(hca,'on'); end
        end
        hold(hca,'off');
        if doTagColor && not(skipTagColor)
          legend(hbox(ihsave),uniquetags{:})
        end
      end
      drawnow
      disp('Linking: XLim, YLim.')
      hlinks = linkprop(ax,{'XLim','YLim'});
      if not(doPlotIntoAxes)
        hlinks.Targets(1).XLim = [xmin xmax] + 5*[-1 1];
        hlinks.Targets(1).YLim = [zmin zmax] + 2*[-1 1];
      end
      
      if nargout == 1
        varargout{1} = ax;
      elseif nargout == 2
        varargout{1} = ax;
        varargout{2} = hbox;
      end
    end
              
    % Data analysis routines, time derivatives, interpolation, etc.
    function out = get_peaks(obj,nPeaks,spacingPeaks,iSpecies,varargin)
      
      doVxlim = 0;
      doVylim = 0;
      doVzlim = 0;
      have_options = 0;
      nargs = numel(varargin);      
      if nargs > 0, have_options = 1; args = varargin(:); end
      while have_options
        l = 1;
        %lower(args{1})
        switch(lower(args{1}))
          case {'vx'}
            doVxlim = 1;
            vxlim = args{2};
            l = 2;         
          case {'vy'}
            doVylim = 1;
            vylim = args{2};
            l = 2;  
          case {'vz'}
            doVzlim = 1;
            vzlim = args{2};
            l = 2;                 
          otherwise
            warning(sprintf('Input ''%s'' not recognized.',args{1}))            
        end        
        args = args(l+1:end);
        if isempty(args), break, end    
      end
      
      nTimes = obj.nt;
      nDists = obj.nd;
      fpeaks = struct;
      for it = 1:nTimes
        for id = 1:nDists{it}
          f = obj.f(it,id,iSpecies);
          ftmp = f.f;
          v = f.v;
          if doVxlim, vxind = union(find(v<vxlim(1)),find(v>vxlim(2)));
          else, vxind = []; end
          if doVylim, vyind = union(find(v<vylim(1)),find(v>vylim(2)));
          else, vyind = []; end
          if doVzlim, vzind = union(find(v<vzlim(1)),find(v>vzlim(2)));
          else, vzind = []; end
          
          %ftmp(vxind,vyind,vzind) = NaN;
          ftmp(vxind,:,:) = NaN;
          ftmp(:,vyind,:) = NaN;
          ftmp(:,:,vzind) = NaN;
          
          for iPeak = 1:nPeaks
            [val,ind] = max(ftmp(:));
           % iPeak
           % val 
           % ind
           
            if val == 0
              continue
            end
            try
            [ix,iy,iz] = ind2sub(size(ftmp),ind);
            ix_ = ix+[-spacingPeaks:spacingPeaks]; ix_(ix_<0) = [];
            iy_ = iy+[-spacingPeaks:spacingPeaks]; iy_(iy_<0) = [];
            iz_ = iz+[-spacingPeaks:spacingPeaks]; iz_(iz_<0) = [];
            ftmp(ix_,iy_,iz_) = NaN;
            catch
              1;
            end
%             if f.v(iz)<0
%               1;
%             end
            fpeaks(iPeak,id,it).vx = f.v(ix);
            fpeaks(iPeak,id,it).vy = f.v(iy);
            fpeaks(iPeak,id,it).vz = f.v(iz);
            fpeaks(iPeak,id,it).x = mean(f.x);
            fpeaks(iPeak,id,it).y = 0;
            fpeaks(iPeak,id,it).z = mean(f.z);
            fpeaks(iPeak,id,it).f = val;
            
            fpeaks(iPeak,id,it).nPeaks = nPeaks;
            fpeaks(iPeak,id,it).spacingPeaks = spacingPeaks;           
            fpeaks(iPeak,id,it).iPeak = iPeak;
            fpeaks(iPeak,id,it).iSpecies = iSpecies;
            fpeaks(iPeak,id,it).timeIteration = obj.iteration(it);
            fpeaks(iPeak,id,it).distId = obj.dists{it}{id};
            fpeaks(iPeak,id,it).dist_x1 = obj.xi1{it}(id);
            fpeaks(iPeak,id,it).dist_x2 = obj.xi2{it}(id);
            fpeaks(iPeak,id,it).dist_z1 = obj.zi1{it}(id);
            fpeaks(iPeak,id,it).dist_z2 = obj.zi2{it}(id);
          end
        end
      end
      out = fpeaks;
    end
    
    % Get simulation meta data and parameters
    function out = get_twpe(obj)
      fileInfo = obj.info_;
      iGroup = find(contains({fileInfo.Groups.Name},'/data'));
      nIter = numel(fileInfo.Groups(iGroup).Groups); % number of iterations
      
      for iIter = 1:nIter
        % /data/00000xxxxx/ 
        % redo to actually find the 
        iAtt = find(contains({fileInfo.Groups(iGroup).Groups(iIter).Attributes.Name},'twpe'));
        time(iIter) = fileInfo.Groups(iGroup).Groups(iIter).Attributes(iAtt).Value;
      end
      out = time;
    end
    function out = get_twci(obj)
      fileInfo = obj.info_;
      iGroup = find(contains({fileInfo.Groups.Name},'/data'));
      nIter = numel(fileInfo.Groups(iGroup).Groups); % number of iterations
      
      for iIter = 1:nIter
        % /data/00000xxxxx/ 
        % redo to actually find the 
        iAtt = find(contains({fileInfo.Groups(iGroup).Groups(iIter).Attributes.Name},'twci'));
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
    function out = get_distlist(obj)
      % needs to be adapted for the species subgroups
      fileInfo = obj.info;
      iterations = obj.iteration;
      % fields structure is the same for all times
      for iIter = 1:numel(iterations)
        dists_iter = cell(numel(fileInfo.Groups(1).Groups(iIter).Groups),1);
        %dists{iIter} = numel(fileInfo.Groups(1).Groups(1).Groups);
        %dists_iter = {fileInfo.Groups(1).Groups(1).Groups.Name};
        split_str = cellfun(@(x) strsplit(x,'/'),{fileInfo.Groups(1).Groups(iIter).Groups.Name},'UniformOutput',false);
        for iout = 1:numel(split_str)
          dists_iter{iout} = split_str{iout}{end};
        end
        dists{iIter} = dists_iter;        
      end
      out = dists;
    end
    function out = get_tags(obj)
      % Get tags, all distributions does not have tags
      fileInfo = obj.info_;
      iGroup = find(contains({fileInfo.Groups.Name},'/data'));
      nIter = numel(fileInfo.Groups(iGroup).Groups); % number of iterations
      all_tags = cell(1,nIter);
      for iIter = 1:nIter
        % /data/00000xxxxx/ 
        % redo to actually find the         
        nDists = numel(fileInfo.Groups(iGroup).Groups(iIter).Groups);
        
        for iDist = 1:nDists
          % Try group attributes first (i saved them a bit differently)
          if isempty(fileInfo.Groups(iGroup).Groups(iIter).Groups(iDist).Attributes)
            all_tags{iIter}{iDist} = '';
          else
            iAtt =   find(contains({fileInfo.Groups(iGroup).Groups(iIter).Groups(iDist).Attributes.Name},'tag')); 
            if isempty(iAtt) % go to subdataset
              iAtt = find(contains({fileInfo.Groups(iGroup).Groups(iIter).Groups(iDist).Datasets(1).Attributes.Name},'tag')); 
              if isempty(iAtt) % no attribute found
                all_tags{iIter}{iDist} = '';
              end
            else
              all_tags{iIter}{iDist} = fileInfo.Groups(iGroup).Groups(iIter).Groups(iDist).Attributes(iAtt).Value;
            end
          end
          if isempty(all_tags{iIter}{iDist})
            all_tags{iIter}{iDist} = ' ';
          end
%           iAtt =      find(contains({fileInfo.Groups(iGroup).Groups(iIter).Groups(iDist).Datasets(1).Attributes.Name},'tag')); 
%           if isempty(iAtt)
%             all_tags{iIter}{iDist} = '';
%           else
%             all_tags{iIter}{iDist} = fileInfo.Groups(iGroup).Groups(iIter).Groups(iDist).Datasets(1).Attributes(iAtt).Value;
%           end
        end
        
      end
      out = all_tags;
      
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
    
    function [x1,x2,z1,z2] = get_locs(obj)
      % Read locations of distributions
      iterations = obj.iteration;
      
      x1 = cell(1,obj.nt);
      x2 = cell(1,obj.nt);
      z1 = cell(1,obj.nt);
      z2 = cell(1,obj.nt);
      
      for iIter = 1:obj.nt
        iter = iterations(iIter);
        str_iter = sprintf('%010.0f',iter);
        % this gives error if there is a jump in dists
        % use obj.dists{iIter}{iDist} instead of num2str(iDist,'%05.0f')?
        for iDist = 1:numel(obj.dists{iIter}) 
          %dataset_name = ['/data/' str_iter '/' num2str(iDist,'%05.0f') '/fxyz'];
          dataset_name = ['/data/' str_iter '/' obj.dists{iIter}{iDist} '/fxyz'];          
          %locs{iIter}{iDist} = [h5readatt(obj.file,dataset_name,'x') h5readatt(obj.file,dataset_name,'z')];          
          try
          xtmp = h5readatt(obj.file,dataset_name,'x')';          
          ztmp = h5readatt(obj.file,dataset_name,'z')';
          catch
           xtmp = [NaN NaN];
           ztmp = [NaN NaN];
          end
          x1{1,iIter}(iDist) = xtmp(1);
          x2{1,iIter}(iDist) = xtmp(2);
          z1{1,iIter}(iDist) = ztmp(1);
          z2{1,iIter}(iDist) = ztmp(2);          
        end
      end      
    end  
    function [v1,v2,nv,dv] = get_v(obj)
      % Read locations of distributions
      iterations = obj.iteration;            
      
      v1 = cell(1,obj.nt);
      v2 = cell(1,obj.nt);
      nv = cell(1,obj.nt);
      dv = cell(1,obj.nt);
      
      for iIter = 1:obj.nt
        iter = iterations(iIter);
        str_iter = sprintf('%010.0f',iter);
        for iDist = 1:numel(obj.dists{iIter}) 
          dataset_name = ['/data/' str_iter '/' obj.dists{iIter}{iDist} '/fxyz'];          
          try
            vtmp = h5readatt(obj.file,dataset_name,'axes')';
          catch
            vtmp = [NaN NaN];
          end
          v1{1,iIter}(iDist,:) = vtmp(:,1);
          v2{1,iIter}(iDist,:) = vtmp(:,end);
          nv{1,iIter}(iDist,:) = size(vtmp,2);
          dv{1,iIter}(iDist,:) = vtmp(:,2)-vtmp(:,1);
        end
      end      
    end  
    function out = nd(obj)
      nd = cellfun(@(x) numel(x),obj.dists,'UniformOutput',false);
      out = nd;
    end
    function out = nx(obj)
      nx = cellfun(@(x) numel(x),obj.dists,'UniformOutput',false);
      out = nx;
    end
    function out = nz(obj)
      nx = cellfun(@(x) numel(x),obj.dists,'UniformOutput',false);
      out = nx;
    end
    function out = nt(obj)
      out = numel(obj.iteration);
    end
    function out = xi(obj)
    end
    function out = zi(obj)
    end
    
    % Phase space distribution
    function out = dataset_str(obj,it,id)
      % PICDist.DATASET_STR
      %   Returns dataset corresponding to obj (object), it (iteration 
      %   index), and id (distribution number)
      %   
      %   dst.dataset_str(1,10)
      %
      %   >> HDF5 dists.h5 
      %   Dataset 'fxyz' 
      %      Size:  101x101x101x6
      %       MaxSize:  101x101x101x6
      %       Datatype:   H5T_IEEE_F64LE (double)
      %       ChunkSize:  []
      %       Filters:  none
      %       FillValue:  0.000000
      %       Attributes:
      %           'x':  182.250000 182.750000 
      %           'z':  2.750000 3.250000 
      %           'ic':  1152065536.000000 1153564672.000000 1150148608.000000 1145405440.000000 0.000000 1125449728.000000 
      %           'vxa':  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 
      %           'vya':  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 
      %           'vza':  0.000000 0.000000 0.000000 0.000000 0.000000 0.000000 
      %           'axes':  101x6 H5T_FLOAT
      
      out = ['/data/' num2str(obj.iteration(it),'%010.0f') '/' num2str(obj.indices{it}(id),'%05.0f') '/fxyz'];
    end
    function varargout = fxyz(obj,it,id,iss,sumdim)
            
      iSpecies = iss;        
      dataset = obj.dataset_str(it,id);
      twci = obj.twci(it);
      twpe = obj.twpe(it);
      info = h5info(obj.file,dataset);      
      nSpecies = info.Dataspace.Size(4);
      datasize = info.Dataspace.Size(1:3);
            
      allAxes = h5readatt(obj.file,dataset,'axes');
      selAxes = allAxes(:,iSpecies);
      nUniqueAxes = size(unique(selAxes','rows'),1);
      comAxes = min(selAxes(:)):min(min(diff(selAxes,1))):max(selAxes(:)); % keep max resolution and max axes value
      newDataSize = [numel(comAxes) numel(comAxes) numel(comAxes)];
      
      dist = zeros(newDataSize);
      %h = setup_subplots(3,2);
      ip = 1;
      for iSpeciesTmp = iSpecies
        % Michael did not divide by the phase space volume, so the units of
        % data_tmp when loaded is #/r3, not #/r3v3, and sum(data_tmp(:)) 
        % gives the density directly
        data_tmp = h5read(obj.file,dataset,[1 1 1,iSpeciesTmp],[datasize,1]);        
        % The phase space volume might be different due to different vaxes
        % or spatial box sizes. Adjust for that here.
        vAxesTmp = allAxes(:,iSpeciesTmp);
        dv = vAxesTmp(2) - vAxesTmp(1);
        
        %hca = h(ip); ip = ip + 1;
        %imagesc(hca,allAxes(:,iSpeciesTmp),allAxes(:,iSpeciesTmp),squeeze(sum(data_tmp)))
        %hca.Title.String = sprintf('sum(data_tmp(:)) = %g',sum(data_tmp(:)));
%         xTmp = h5readatt(obj.file,dataset,'x')';
%         zTmp = h5readatt(obj.file,dataset,'z')';
%         dx = diff(xTmp);
%         dz = diff(zTmp);
%         volPhaseSpace = dv.^3*dx*dz;
%dv = 1;
        volVelocitySpace = dv.^3;
        data_tmp = data_tmp/volVelocitySpace; % #/r3v3
        
        if nUniqueAxes > 1 % Resample/interpolate to common axes.
          [Vx,Vy,Vz] = meshgrid(allAxes(:,iSpeciesTmp),allAxes(:,iSpeciesTmp),allAxes(:,iSpeciesTmp));
          [Vxq,Vyq,Vzq] = meshgrid(comAxes,comAxes,comAxes);
          data_tmp = interp3(Vx,Vy,Vz,data_tmp,Vxq,Vyq,Vzq);
          data_tmp(isnan(data_tmp)) = 0;
          new_dv = comAxes(2)-comAxes(1);
          %new_dv = 1; % ???
          data_tmp = data_tmp*new_dv.^3; % dn
        end       
        %hca = h(ip); ip = ip + 1;
        %imagesc(hca,comAxes,comAxes,squeeze(sum(data_tmp)))
        %hca.Title.String = sprintf('sum(data_tmp(:)) = %g',sum(data_tmp(:)));
        dist = dist + data_tmp;
        
      end
%      hlinks = linkprop(h,{'XLim','YLim','CLim'});
      
      % Reduce distribution
      if exist('sumdim','var') && any(sumdim == [1 2 3])
        for isum = numel(sumdim):-1:1
          dist = squeeze(sum(dist,sumdim(isum)))*(comAxes(2)-comAxes(1)); % fdv
        end
      end
      
      x = h5readatt(obj.file,dataset,'x');
      z = h5readatt(obj.file,dataset,'z');
      
      if nargout == 1        
        out.twci = twci;
        out.twpe = twpe;
        out.f = dist;
        out.v = comAxes;
        out.dv = comAxes(2)-comAxes(1);
        out.x = x;
        out.z = z;
        varargout{1} = out;
      elseif nargout == 3
        varargout{1} = comAxes;
        varargout{2} = comAxes;
        varargout{1} = dist;
      end
    end
    function varargout = f(obj,it,id,iss,varargin)
            
      doPar = 0;
      have_options = 0;
      if not(isempty(varargin))
        args = varargin;
        have_options = 1;
      end
      while have_options
        l = 1;
        switch lower(args{1})
          case 'par'
            doPar = 1;
            l = 1;
        end
        args = args((1+l):end);
      end
      
      iSpecies = iss;        
      dataset = obj.dataset_str(it,id);
      twci = obj.twci(it);
      twpe = obj.twpe(it);
      info = h5info(obj.file,dataset);      
      nSpecies = info.Dataspace.Size(4);
      datasize = info.Dataspace.Size(1:3);
            
      allAxes = h5readatt(obj.file,dataset,'axes');
      selAxes = allAxes(:,iSpecies);
      nUniqueAxes = size(unique(selAxes','rows'),1);
      comAxes = min(selAxes(:)):min(min(diff(selAxes,1))):max(selAxes(:)); % keep max resolution AND max axes value (this can increase the number of bins)
      newDataSize = [numel(comAxes) numel(comAxes) numel(comAxes)];
      
      dist = zeros(newDataSize);
      for iSpeciesTmp = iSpecies
        data_tmp = h5read(obj.file,dataset,[1 1 1,iSpeciesTmp],[datasize,1]);
        % The phase space volume might be different due to different vaxes
        % or spatial box sizes. Adjust for that here. W/hen summing up the
        % untreated distribution, one gets the density, which means the
        % data loaded is actually already multiplied with the velocity volume
        % element, but not the space volume element. Need to undo that here
        % so that distributions with different velocity axes can be added 
        % together.
        vAxesTmp = allAxes(:,iSpecies);
        dv = vAxesTmp(2) - vAxesTmp(1);
        %xTmp = h5readatt(obj.file,dataset,'x')';
        %zTmp = h5readatt(obj.file,dataset,'z')';
        %dx = diff(xTmp);
        %dz = diff(zTmp);
        volVelocitySpace = dv.^3;
        data_tmp = data_tmp/volVelocitySpace;
        
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
          dist_sum{sumdim} = squeeze(sum(dist,sumdim(isum)))*(comAxes(2)-comAxes(1)); % fdv
        end        
      end      
      fx = reshape(squeeze(sum(sum(dist,3),2))*(comAxes(2)-comAxes(1))^2,[size(dist,1),1]); % fdv
      fy = reshape(squeeze(sum(sum(dist,3),1))*(comAxes(2)-comAxes(1))^2,[size(dist,1),1]); % fdv
      fz = reshape(squeeze(sum(sum(dist,2),1))*(comAxes(2)-comAxes(1))^2,[size(dist,1),1]); % fdv
      
      x = h5readatt(obj.file,dataset,'x');
      z = h5readatt(obj.file,dataset,'z');
      
      if doPar
        vpar = comAxes;        
      end
      
      out.twci = twci;
      out.twpe = twpe;
      out.f = dist;
      out.fyz = dist_sum{1};
      out.fxz = dist_sum{2};
      out.fxy = dist_sum{3};      
      out.fx = fx;
      out.fy = fy;
      out.fz = fz;
      out.v = comAxes;
      out.dv = comAxes(2)-comAxes(1);
      out.x = x;
      out.z = z;
      varargout{1} = out;
    end
    function varargout = fs(obj,it,id,iss,comp,varargin)
      iSpecies = iss; 
      if ischar(id) && strcmp(id,':')
        id = 1:obj.nd{1};
      end
      
      for iDist = 1:numel(id)
        iSpeciesTmp = iSpecies;
        ff = fxyz(obj,it,iDist,iSpeciesTmp);
               
        switch comp
          case 'x'
            fred = sum(sum(ff.f,3),2)*ff.dv^2;
          case 'y'
            fred = sum(sum(ff.f,3),1)*ff.dv^2;
          case 'z'
            fred = sum(sum(ff.f,2),1)*ff.dv^2;
        end
        if iDist == 1 % initialize arrays
          dist = zeros(numel(fred),numel(id));
        end
        dist(:,iDist) = fred;
        x(:,iDist) = sum(ff.x)/2;
        z(:,iDist) = sum(ff.z)/2;        
        x1(:,iDist) = ff.x(1);
        z1(:,iDist) = ff.z(1);
        x2(:,iDist) = ff.x(2);
        z2(:,iDist) = ff.z(2);
      end
      % only makes sense if distributions are connected
      arclength = [0 cumsum(sqrt(diff(x).^2+diff(z).^2))];
      arc01 = arclength(find(abs(z)==min(abs(z))));
      arclength_z0 = arclength - arc01(1);

      out.x = x; 
      out.z = z; 
      out.x1 = x1; 
      out.x2 = x2; 
      out.z1 = z1; 
      out.z2 = z2; 
      out.v = ff.v;
      out.f = dist;
      out.arc = arclength;      
      out.arc_z0 = arclength_z0;
      varargout{1} = out;
    end
    function out = fx(obj,it,id,iss,varargin)
      out = fred(obj,it,id,iss,'x',varargin);
    end
    function out = fy(obj,it,id,iss,varargin)
      out = fred(obj,it,id,iss,'y',varargin);
    end
    function out = fz(obj,it,id,iss,varargin)
      out = fred(obj,it,id,iss,'z',varargin);
    end
    function varargout = fpar(obj,it,id,iss,varargin)
      iSpecies = iss;
      %twci = obj.twci(it);
      %twpe = obj.twpe(it);  
      if ischar(id) && strcmp(id,':')
        id = 1:obj.nd{1};
      end
      dataset1 = ['/data/' num2str(obj.iteration(it(1)),'%010.0f') '/' num2str(obj.indices{it(1)}(id(1)),'%05.0f') '/fpar'];      
      allAxes = h5readatt(obj.file,dataset1,'axes');
      datasize = size(allAxes);
      v = allAxes(:,iSpecies);
      dist = zeros(numel(v),numel(id));
      %x = ;
      %z = ;
      for iDist = 1:numel(id)
        iSpeciesTmp = iSpecies;
        dataset = ['/data/' num2str(obj.iteration(it(1)),'%010.0f') '/' num2str(obj.indices{it(1)}(id(iDist)),'%05.0f') '/fpar'];        
        data_tmp = h5read(obj.file,dataset,[1,iSpeciesTmp],[datasize(1),1]);
        x_tmp = h5readatt(obj.file,dataset,'x');
        z_tmp = h5readatt(obj.file,dataset,'z');
        
        dist(:,iDist) = data_tmp;
        x(:,iDist) = sum(x_tmp)/2;
        z(:,iDist) = sum(z_tmp)/2;        
        x1(:,iDist) = x_tmp(1);
        z1(:,iDist) = z_tmp(1);
        x2(:,iDist) = x_tmp(2);
        z2(:,iDist) = z_tmp(2);
      end
      % only makes sense if distributions are connected
      arclength = [0 cumsum(sqrt(diff(x).^2+diff(z).^2))];
      arc01 = arclength(find(abs(z)==min(abs(z))));
      arclength_z0 = arclength - arc01(1);

      out.x = x; 
      out.z = z; 
      out.x1 = x1; 
      out.x2 = x2; 
      out.z1 = z1; 
      out.z2 = z2; 
      out.v = v;
      out.f = dist;
      out.arc = arclength;      
      out.arc_z0 = arclength_z0;
      varargout{1} = out;
    end
    function varargout = fred(obj,it,id,iss,comp,varargin)
      iSpecies = iss;
      %twci = obj.twci(it);
      %twpe = obj.twpe(it);  
      if ischar(id) && strcmp(id,':')
        id = 1:obj.nd{1};
      end
      dataset1 = ['/data/' num2str(obj.iteration(it(1)),'%010.0f') '/' num2str(obj.indices{it(1)}(id(1)),'%05.0f') '/f' comp];      
      allAxes = h5readatt(obj.file,dataset1,'axes');
      datasize = size(allAxes);
      v = allAxes(:,iSpecies);
      dist = zeros(numel(v),numel(id));
      %x = ;
      %z = ;
      for iDist = 1:numel(id)
        iSpeciesTmp = iSpecies;
        dataset_base = ['/data/' num2str(obj.iteration(it(1)),'%010.0f') '/' num2str(obj.indices{it(1)}(id(iDist)),'%05.0f')];        
        dataset = ['/data/' num2str(obj.iteration(it(1)),'%010.0f') '/' num2str(obj.indices{it(1)}(id(iDist)),'%05.0f') '/f' comp];        
        data_tmp = h5read(obj.file,dataset,[1,iSpeciesTmp],[datasize(1),1]);
        %if strcmp(comp,'par')
          %x_tmp = h5readatt(obj.file,dataset,'x');
          %z_tmp = h5readatt(obj.file,dataset,'z');
          x_tmp = h5readatt(obj.file,[dataset_base '/fxyz'],'x');
          z_tmp = h5readatt(obj.file,[dataset_base '/fxyz'],'z');
        %else          
        %end        
        B_tmp = h5readatt(obj.file,dataset_base,'B');
        E_tmp = h5readatt(obj.file,dataset_base,'E');
        
          
        dist(:,iDist) = data_tmp;
        x(:,iDist) = sum(x_tmp)/2;
        z(:,iDist) = sum(z_tmp)/2;        
        x1(:,iDist) = x_tmp(1);
        z1(:,iDist) = z_tmp(1);
        x2(:,iDist) = x_tmp(2);
        z2(:,iDist) = z_tmp(2);
        Bx(:,iDist) = B_tmp(1);
        By(:,iDist) = B_tmp(2);
        Bz(:,iDist) = B_tmp(3);
        Ex(:,iDist) = E_tmp(1);
        Ey(:,iDist) = E_tmp(2);
        Ez(:,iDist) = E_tmp(3);
      end
      % only makes sense if distributions are connected
      arclength = [0 cumsum(sqrt(diff(x).^2+diff(z).^2))];
      % does this needs to be interpolated?
      if 1 % closest
        arc01 = arclength(find(abs(z)==min(abs(z))));
      else % interpolation
        i1 = find(z>0,1,'last');
        i2 = find(z<0,1,'first');
        arc01 = interp1(z([i1 i2]),arclength([i1 i2]),0);
      end
      %arc01 = arclength(find(abs(z)==min(abs(z))));      
      %z0_offs = z(find(abs(z)==min(abs(z))));
      arclength_z0 = arclength - arc01(1);

      out.x = x; 
      out.z = z; 
      out.x1 = x1; 
      out.x2 = x2; 
      out.z1 = z1; 
      out.z2 = z2; 
      out.v = v;
      out.f = dist;
      out.arc = arclength;      
      out.arc_z0 = arclength_z0;
      out.Bx = Bx;
      out.By = By;
      out.Bz = Bz;
      out.Ex = Ex;
      out.Ey = Ey;
      out.Ez = Ez;
      
      varargout{1} = out;
    end
    
    function varargout = rebin(obj,it,id,iss,varargin)
      % PICDist.rebin Rebin distribution.
      %
      % Options
      %   space - logarithmically spaced energies with deltaE/E constant
      %   vabs - absolute value of v
      %   E - 2*v^2/m
      %   proj - 1D along a line
      %   cart - 3D
      %   sphere - spherical            
      
      have_options = 0;
      if not(isempty(varargin))
        args = varargin;
        have_options = 1;
      end
      while have_options
        l = 1;
        switch lower(args{1})
          case 'proj'
            doProj = 1;
            l = 1;
          case 'vabs'
            rebin_option = 'vabs';
            vgrid = varargin{2};
            l = 2;
          case 'space'
            rebin_option = 'space';
            l = 1;
            
          case 'sphere'
            
          case 'cart'
        end
        args = args((1+l):end);
        if isempty(args); break; end
      end
      
      % just assume for now every other is proton and electrons, needs to
      % be included in obj.properties
      masses = [25 1 25 1 25 1 25 1 25 1]/25; 
      mass = masses(iss);
      % 
      id_count = 0;
      for id = 1:numel(ids)
        id_count = id_count + 1;
      
        f = obj.fxyz(it,id,iss);
        [VX,VY,VZ] = ndgrid(f.v,f.v,f.v); % 3D mesh of grid
        VABS = sqrt(VX.^2+VY.^2+VZ.^2); % magnitude of velocity or all cells
        ENERGY = mass*VABS.^2/2;
        d3v = (f.v(2)-f.v(1))^3; % phase space volume of each bin, so to 
                                 % get the density of each bin do f*d3v,
                                 % and to get the total density, do
                                 % sum(f(:)*d3v)
        
        switch rebin_option
          case 'vabs'
            if isempty(vgrid) % if not supplied, define here based on max(VABS(:))
              nv = 70;
              vgrid = tocolumn(linspace(0,max(VABS(:)),nv));
            end
            dv = vgrid(2:nv) - vgrid(1:(nv-1));
            vcenter = vgrid(1:(nv-1)) + dv;
            
            % For givin output in multiple units
            Egrid = mass*vgrid.^2/2;
            %Ecenter = mass*vcenter.^2/2;
            Eminus = mass*vgrid(1:(nv-1)).^2/2;
            Eplus = mass*vgrid(2:(nv)).^2/2;
            Ecenter = (Eplus + Eminus)/2;
            Edelta = Eplus - Eminus;
            
            [count edges mid loc] = histcn(VABS(:),vgrid);
            volume_of_shell = count*d3v;
            % sum(fd3v)
            % units: [f][v3] = l-3 v-3 v3 = l-3 (# per volume)
            % So this is number of particles per volume = density
            [accum_f edges mid loc] = histcn(VABS(:),vgrid,'AccumData',f.f(:)*d3v);
            % sum(fvd3v)
            % units [f][v][v3] = l-3 v-3 v v3 = vl-3 = l s-1 l-1 = s-1 l-2 (# per area and time)
            [accum_fv edges mid loc] = histcn(VABS(:),vgrid,'AccumData',f.f(:).*VABS(:)*d3v);            
            % sum(fvEd3v)
            % units [f][v][E][v3] = l-3 v-3 v eV v3 = v eV l-3 = l s-1 eV l-1 = eV s-1 l-2 (energy per area and time)
            [accum_fvE edges mid loc] = histcn(VABS(:),vgrid,'AccumData',f.f(:).*VABS(:).*ENERGY(:)*d3v);
            % Since we used equally spaced v's, the new phase space volumes
            % (which are shells), become increasingly big further out
            % (larger vabs): d3v = v2*dv*dOmega. Omega is the solid angle
            % and is the same (4*pi for a sphere) for all new bins since we 
            % only binned by |v|.z.
            fout.volume = volume_of_shell;            
            fout.vedges = vgrid;
            fout.dv = dv;
            fout.Eedges = Egrid;
            fout.dE = Edelta;
            fout.vcenter = vcenter;
            fout.Ecenter = Ecenter;
            fout.accum_fd3v = accum_f;
            % units 
            % [solid angle] = sr (steradians)
            % [dv] = l s-1
            % [E] = eV
            fout.dfdv = accum_f/(4*pi)./dv; % m-3/(m/s)= s/m4
            fout.dfdE = accum_f/(4*pi)./Edelta; % s/(m^4 sr eV)
            fout.dfvdv = accum_fv/(4*pi)./dv; % (s-1 m-2)/(m/s) = 1/m3
            fout.dfvdE = accum_fv/(4*pi)./Edelta; % 1/(sr eV) = 1/(s m2 sr eV)
            fout.dfvEdv = accum_fvE/(4*pi)./dv; % 1/(m/s)
            fout.dfvEdE = accum_fvE/(4*pi)./Edelta; % (eV s-1 m-2)/(sr eV) = eV/(s m2 sr eV)
          case 'proj'
            
          case 'space'
            
        end
      
      end
      varargout{1} = fout;
    end
    % dEF - differential energy flux
    function varargout = dEF(obj,it,ids,iss,nE)
      % out = dEF(obj,it,id,iss,nE)
      
      % set up 
      %energy_edges = logspace(-3,log10(1.0*max(f.v.^2/2)),30);
      nEnergy = nE;
      %energy_edges = logspace(-3,log10(1*max(f.v.^2/2)),nEnergy);
      energy_edges = logspace(-3,1,nEnergy);
      %energy_edges = linspace(1e-4,1e1,nEnergy);
      
      % initialize variable(s)
      dEF_all = zeros(numel(ids),nEnergy-1);
      
      %  loop through ids
      id_count = 0;
      for id = 1:numel(ids)
        id_count = id_count + 1;
      
        f = obj.fxyz(it,id,iss);
        
        
        %energy = logspace(-3,log10(max(f.v.^2)),nEnergy);
        [VX,VY,VZ] = ndgrid(f.v,f.v,f.v);
        VV2 = VX.^2 + VY.^2 + VZ.^2;
        fvv = f.f.*VV2;
        ENERGY = VV2/2;        
        [N,EDGES,BIN] = histcounts(ENERGY(:),energy_edges);
        f_energy_edges = tocolumn(EDGES);
        f_denergy = diff(f_energy_edges);
        f_energy_centers = tocolumn((EDGES(2:end)+EDGES(1:end-1))*0.5);
        nbins = (numel(energy_edges)-1);
        for ibin = 1:nbins
          ind_bin = find(BIN==ibin);      
          f_dist_tmp = f.f(ind_bin).*VV2(ind_bin).^2/2; % IS THIS RIGHT NOW ???
          f_dist_mean(ibin,1) = mean(f_dist_tmp);
          %f_dist_sum(ibin,1) = sum(f_dist_tmp)/f_denergy(ibin);
        end

  
        dEF_all(id,:) = f_dist_mean; % NEED TO MULTIPLY WITH V^4/2
        x_all(id) = mean(f.x);
        z_all(id) = mean(f.z);
      end
      if nargout == 1
        varargout{1}.energy = f_energy_centers;
        varargout{1}.energy_edges = f_energy_edges;
        varargout{1}.dEF = dEF_all;
        varargout{1}.x = x_all;
        varargout{1}.z = z_all;
      elseif nargout == 2        
        varargout{2}.energy_edges = f_energy_edges;
        varargout{3}.dEF = dEF_all;
      elseif nargout == 3
        varargout{1}.energy = f_energy_centers;
        varargout{2}.energy_edges = f_energy_edges;
        varargout{3}.dEF = dEF_all;
      end
    end
    function varargout = omni(obj,it,ids,iss,nE)
      % out = dEF(obj,it,id,iss,nE)
      
      % set up 
      %energy_edges = logspace(-3,log10(1.0*max(f.v.^2/2)),30);
      nEnergy = nE;
      %energy_edges = logspace(-3,log10(1*max(f.v.^2/2)),nEnergy);
      energy_edges = logspace(-3,5,nEnergy);
      denergy = diff(energy_edges);
      %energy_edges = linspace(1e-4,1e1,nEnergy);
      
      % Initialize variable(s)
      f_all = zeros(numel(ids),nEnergy-1);
      dPF_all = zeros(numel(ids),nEnergy-1);
      dEF_all = zeros(numel(ids),nEnergy-1);
      
      
      %  loop through ids
      id_count = 0;
      for id = 1:numel(ids)
        id_count = id_count + 1;
      
        f = obj.fxyz(it,id,iss);
        
        
        %energy = logspace(-3,log10(max(f.v.^2)),nEnergy);
        [VX,VY,VZ] = ndgrid(f.v,f.v,f.v);
        VV2 = VX.^2 + VY.^2 + VZ.^2;
        fvv = f.f.*VV2;
        ENERGY = VV2/2;        
        [N,EDGES,BIN] = histcounts(ENERGY(:),energy_edges);
        f_energy_edges = tocolumn(EDGES);
        f_denergy = diff(f_energy_edges);
        f_energy_centers = tocolumn((EDGES(2:end)+EDGES(1:end-1))*0.5);
        nbins = (numel(energy_edges)-1);
        for ibin = 1:nbins
          ind_bin = find(BIN==ibin);      
          E_dist_tmp = f.f(ind_bin).*VV2(ind_bin).^2/2/denergy(ibin); % IS THIS RIGHT NOW ??? No, need to divide by bin size i guess
          E_dist_mean(ibin,1) = mean(E_dist_tmp);
          P_dist_tmp = f.f(ind_bin).*VV2(ind_bin)/denergy(ibin); % IS THIS RIGHT NOW ???
          P_dist_mean(ibin,1) = mean(P_dist_tmp);
          f_dist_tmp = f.f(ind_bin); % IS THIS RIGHT NOW ???
          f_dist_mean(ibin,1) = mean(f_dist_tmp);
          %f_dist_sum(ibin,1) = sum(f_dist_tmp)/f_denergy(ibin);
        end

  
        dEF_all(id,:) = E_dist_mean; % NEED TO MULTIPLY WITH V^4/2
        dPF_all(id,:) = P_dist_mean; % NEED TO MULTIPLY WITH V^4/2
        f_all(id,:) = f_dist_mean; % NEED TO MULTIPLY WITH V^4/2
        x_all(id) = mean(f.x);
        z_all(id) = mean(f.z);
      end
      if nargout == 1
        varargout{1}.energy = f_energy_centers;
        varargout{1}.energy_edges = f_energy_edges;
        varargout{1}.f = f_all;
        varargout{1}.dEF = dEF_all;
        varargout{1}.dPF = dPF_all;
        varargout{1}.x = x_all;
        varargout{1}.z = z_all;
      elseif nargout == 2        
        varargout{2}.energy_edges = f_energy_edges;
        varargout{3}.dEF = dEF_all;
      elseif nargout == 3
        varargout{1}.energy = f_energy_centers;
        varargout{2}.energy_edges = f_energy_edges;
        varargout{3}.dEF = dEF_all;
      end
    end    
    
    function varargout = reduce_1d(obj,depvar,x0,z0,vaxes,iSpecies,varargin)
      % reduce_1d(obj,depvar,x0,z0,vaxes,iSpecies) Make f(x,vx),f(z,vx), etc...
      %
      % apply time indices and xlim zlim outside before
      % do all vdims, doesnt take much more time anyways
      % return as structure array, arr(t,x_or_z)
      %
      %   depvar - 'x' or 'z', variable that is plotted on x-axis, for 
      %     example if depvar = 'x', the output array has dimensions 
      %     (nt,nz), and the matrices has dimensions (nv,nx):
      %     arr(1,1) = f(t=t1,x=xvals,z=zvals(1))      
      %   x0 - needs to match exact distribution centers
      %   z0 - needs to match exact distribution centers
      %   vaxes - common velocity axes to interpolate fields to, for
      %     plotting, if vaxes is empty, just use the axes of the first
      %     distribution
      
      method = 'exact'; % only keep distributions with exact center coordinates
      doV2 = 0;
      doPar = 0;
      
      % Check for additional input
      args = varargin;
      nargs = numel(args);
      have_options = nargs > 0;
      while have_options
        switch(lower(args{1}))
          case {'v2','vabs'} % calculate v^2 (similar to dEF), do not include by default because it might take extra time
            l = 1;            
            doV2 = 1;
          case {'vpar','par'}
            l = 2;
            doPar = 1;
            B = args{2};
          otherwise
            error(sprintf('Can not recognize flag ''%s'' ',args{1}))
        end
        args = args((l+1):end);
        if isempty(args), break, end
      end
      
      % First check all pairs of x and z, if one of them is unique, for
      % example the z = 0 row, then make 2D matrix, that is prepared for
      % plotting.
      
      %ds = obj.xfind(x0).zfind(z0); % moved this below
      
      
      nx = numel(x0);
      nz = numel(z0);
      nv = numel(vaxes);
      
      if strcmp(depvar,'x')
        ns = nz;
      else
        ns = nx;
      end
      
      % for vabs distribution
      vaxes_pos = vaxes(vaxes>=0);
      nv_pos = numel(vaxes_pos-1); % vaxes_pos are the edges of the bins, so the number of bins is one less
      
      % save results in structure array, one structure for each time
      % also save time, but for now only the iteration is saved
      fout = struct('time',{},'iter',{},'x',{},'z',{},'v',{},'fvx',{},'fvy',{},'fvz',{});
         
      for itime = 1:obj.nt % loop through times
        dst = obj.twcifind(obj.twci(itime));
        if isempty(dst.xi1), continue; end
        
        for is = 1:ns % loop through x or z
          if strcmp(depvar,'x')
            ds = dst.xfind(x0).zfind(z0(is)); % loop through zvals
          else
            ds = dst.xfind(x0(is)).zfind(z0); % loop through xvals
          end
          if isempty(ds.xi1), continue; end
        
%         t_count = 0; 
          t_count = itime;
          tic
          disp(sprintf('it = %g/%g, is = %g/%g',itime,obj.nt,is,ns))          
          ids = ds.indices{1}; % this produces error if empty, fix this, because it might interrupt a long run
          
          nd = numel(ids);
          if isempty(ids)
            continue
          %else
          %  t_count = t_count + 1; % this not necessary because empty
          %  entries in the structure matrix is ok
          end
          % f(x,z,vx,vy,vz), typically z or x is single value, in future also
          % make possible to make along a line, for example separatrix
          f_vx_arr = zeros(nd,nv);
          f_vy_arr = zeros(nd,nv);
          f_vz_arr = zeros(nd,nv);
          f_vabs_sum_arr = zeros(nd,nv_pos);
          f_vabs_mean_arr = zeros(nd,nv_pos);
          x_arr = zeros(nd,1);
          z_arr = zeros(nd,1);

          % implement later:
          % sort indices, in case they wouldnt already be. otherwise plotting
          % might be strange if one doesnt specify x or z


          id_count = 0;
          for id = 1:numel(ids)
            id_count = id_count + 1;

            f_tmp = ds.f(1,id,iSpecies);

            % Never mind entering the vaxes, just take whatever comes out
            % here, which depends on iSpecies. But what if there are
            % different vaxes for different disitributions, times? Can one go
            % trough and check that first?          

            x_arr(id) = mean(f_tmp.x);
            z_arr(id) = mean(f_tmp.z);          
            f_vx_arr(id,:) = interp1(f_tmp.v,sum(f_tmp.fxy,2),vaxes);
            f_vy_arr(id,:) = interp1(f_tmp.v,sum(f_tmp.fxy,1),vaxes);
            f_vz_arr(id,:) = interp1(f_tmp.v,sum(f_tmp.fxz,1),vaxes);

            if doV2
              [VX,VY,VZ] = ndgrid(f_tmp.v,f_tmp.v,f_tmp.v);
              VABS = sqrt(VX.^2 + (VY).^2 + VZ.^2);
              %fvv = f_tmp.f.*VV2; % this is not used..?
              %ENERGY = VV2/2;
              [N,EDGES,BIN] = histcounts(VABS(:),vaxes_pos);
              f_vabs_edges = tocolumn(EDGES);
              f_dvabs = diff(f_vabs_edges);
              f_vabs_centers = tocolumn((EDGES(2:end)+EDGES(1:end-1))*0.5);
              nbins = nv_pos;
              for ibin = 1:nbins
                %try
                %ibin
                ind_bin = find(BIN==ibin);      
                f_dist_tmp = f_tmp.f(ind_bin);
                f_mean_tmp = nanmean(f_dist_tmp);
                f_sum_tmp = nansum(f_dist_tmp);
                if f_mean_tmp>0
                  1;
                end
                if not(isempty(f_mean_tmp))
                  f_vabs_mean_arr(id,ibin) = f_mean_tmp;
                end
                if not(isempty(f_sum_tmp))
                  f_vabs_sum_arr(id,ibin) = f_sum_tmp;%/f_dvabs(ibin);
                end
                %catch
                  1;
                %end
              end
            end    
            if doPar
              [VX,VY,VZ] = ndgrid(f_tmp.v,f_tmp.v,f_tmp.v);
            end
          end
          [x_arr_sorted,i_sorted] = sort(x_arr);
          
          fout(t_count,is).time = ds.twci(1);
          fout(t_count,is).iter = ds.iteration(1);
          fout(t_count,is).x = x_arr(i_sorted);
          fout(t_count,is).z = z_arr;
          fout(t_count,is).v = vaxes;
          fout(t_count,is).fvx = f_vx_arr(i_sorted,:);
          fout(t_count,is).fvy = f_vy_arr(i_sorted,:);
          fout(t_count,is).fvz = f_vz_arr(i_sorted,:);
          fout(t_count,is).vabs_edges = vaxes_pos;
          [X,V] = meshgrid(fout(t_count,is).x,vaxes_pos);
          fout(t_count,is).vabs_edges = vaxes_pos;
          fout(t_count,is).vabs_mat = V; % multiply f_vabs_mean_arr(i_sorted,:) with vabs_mat.^4/2 to get DEF        
          fout(t_count,is).fvabssum = f_vabs_sum_arr(i_sorted,:);
          fout(t_count,is).def = f_vabs_mean_arr(i_sorted,:).*V'.^4/2;

          fout(t_count,is).fvabsmean = f_vabs_mean_arr(i_sorted,:);

          toc
        end
      end
      varargout{1} = fout;
    end
    function varargout = reduce_1d_new(obj,depvar,iSpecies,vaxes,varargin)
      % reduce_1d(obj,depvar,vaxes,iSpecies) Make f(x,vx),f(z,vx), etc...
      %
      % apply time indices and xlim zlim outside before
      % do all vdims, doesnt take much more time anyways
      % return as structure array, arr(t,x_or_z)
      %
      %   depvar - 'x' or 'z', variable that is plotted on x-axis, for 
      %     example if depvar = 'x', the output array has dimensions 
      %     (nt,nz), and the matrices has dimensions (nv,nx):
      %     arr(1,1) = f(t=t1,x=xvals,z=zvals(1))      
      %   vaxes - common velocity axes to interpolate fields to, for
      %     plotting, if vaxes is empty, just use the axes of the first
      %     distribution
            
      doV2 = 0;
      doE = 0;
      doPar = 0;
      doVabs = 0;
      doPitch = 0;
      doInterp = 0;
      
      % Check for additional input
      args = varargin;
      nargs = numel(args);
      have_options = nargs > 0;
      while have_options
        switch(lower(args{1}))
          case {'v2','vabs'} % calculate v^2 (similar to dEF), do not include by default because it might take extra time
            l = 1;            
            doV2 = 1;
          case {'vpar','par'}
            l = 2;
            doPar = 1;
            B = args{2}; % {Bx,By,Bz} 
          case {'pitch'}
            l = 2;
            doPitch = 1;
            B = args{2};
          case 'interp'            
            nInterp = args{2};
            if nInterp == 0
              doInterp = 0;
            else
              doInterp = 1;              
            end
            l = 2;
          otherwise
            error(sprintf('Can not recognize flag ''%s'' ',args{1}))
        end
        args = args((l+1):end);
        if isempty(args), break, end
      end
      
      % First check all pairs of x and z, if one of them is unique, for
      % example the z = 0 row, then make 2D matrix, that is prepared for
      % plotting.
      
      %ds = obj.xfind(x0).zfind(z0); % moved this below
      
      if isempty(vaxes)
        [v1,v2,nv,dv] = obj.get_v;
        v1 = v1{1}(:,iSpecies);
        v2 = v2{1}(:,iSpecies);
        nv = nv{1}(:);
        dv = dv{1}(:,iSpecies);
        vaxes = min(v1(:)):min(dv(:)):max(v2(:));
      end
            
      %nd = sum([obj.nd{:}]);
      nv = numel(vaxes);
      
      n_vabs = 50;
      n_vpar = 101;
      n_pitch = 18;
      n_v_pitch = 20;
      n_E_pitch = 20;
      % save results in structure array, one structure for each time
      % also save time, but for now only the iteration is saved
      %fout = struct('time',{},'iter',{},'x',{},'z',{},'v',{},'fvx',{},'fvy',{},'fvz',{});
                           
      for itime = 1:obj.nt % loop through times
        ds = obj.twcifind(obj.twci(itime));
        if isempty(ds.xi1), continue; end       
        ids = ds.indices{1};
        nd = obj.nd{1};
        f_vx_arr = zeros(nd,nv);
        f_vy_arr = zeros(nd,nv);
        f_vz_arr = zeros(nd,nv);
        f_vabs_arr = zeros(nd,n_vabs);
        f_vpar_arr = zeros(nd,n_vpar);
        f_pitch_arr = zeros(nd,n_v_pitch,n_pitch);
        f_pitchE_arr = zeros(nd,n_E_pitch,n_pitch);
        def_E_arr = zeros(nd,n_vabs);
        dpf_E_arr = zeros(nd,n_vabs);
        t_arr = zeros(nd,1);
        x_arr = zeros(nd,1);
        z_arr = zeros(nd,1);
          
        id_count = 0;        
        fprintf('id = %4.0f/%4.0f\n',0,numel(ids)) % display progress              
        for id = 1:numel(ids) % loop through distributions
          if mod(id,1) == 0, fprintf([repmat('\b', 1, 10) '%4.0f/%4.0f\n'],id,numel(ids)); end % display progress
          %disp(sprintf('it = %g/%g, id = %g/%g',itime,obj.nt,id,numel(ids)))
%         t_count = 0; 
          t_count = itime;
          %tic           
          
          id_count = id_count + 1;
          f = ds.f(1,id,iSpecies);       

          t_arr(id_count) = ds.twci;
          x_arr(id_count) = mean(f.x);
          z_arr(id_count) = mean(f.z);          
          dv_old = f.v(2)-f.v(1);
          dv_new = vaxes(2)-vaxes(1);
          % the interpolation may change the summed up (phase space)
          % density? test:
          if 0 
            %%
            nend = 100;
            aa = 0:nend;
            ff = 0:100;
            aa_new = 0:7:nend;
            ff_new = interp1(aa,ff,aa_new);
            sum(ff)/(sum(ff_new)*(aa_new(2)-aa_new(1))/(aa(2)-aa(1)))
          end
          f_vx_arr(id_count,:) = interp1(f.v,sum(f.fxy,2)*dv_old,vaxes);%*dv_new/dv_old;
          f_vy_arr(id_count,:) = interp1(f.v,sum(f.fxy,1)*dv_old,vaxes);%*dv_new/dv_old;
          f_vz_arr(id_count,:) = interp1(f.v,sum(f.fxz,1)*dv_old,vaxes);%*dv_new/dv_old;
          %f_vabs_arr(id_count,:) = interp1(f_tmp.v,sum(f_tmp.fxz,1),vaxes);
          


          % Just do them together
          if any([doVabs doE])
            % for vabs distribution
            dv = f.v(2)-f.v(1);
            d3v = dv.^3;
            [VX,VY,VZ] = ndgrid(f.v,f.v,f.v);
            VABS = sqrt(VX.^2 + (VY).^2 + VZ.^2);
            v_abs_grid = linspace(0,max(VABS(:)),n_vabs+1);
            dv_abs = v_abs_grid(2)-v_abs_grid(1);
            v_abs_center = v_abs_grid(1:(end-1)) + dv_abs;
            
            masses = [25 1 25 1 25 1];
            mass = masses(iSpecies(1))/masses(1);
            ENERGY = mass*VABS.^2/2;
            
            % For givin output in multiple units
            Egrid = mass*v_abs_grid.^2/2;
            %Ecenter = mass*vcenter.^2/2;
            Eminus = mass*v_abs_grid(1:(n_vabs)).^2/2;
            Eplus = mass*v_abs_grid(2:(n_vabs+1)).^2/2;
            Ecenter = (Eplus + Eminus)/2;
            Edelta = Eplus - Eminus;
            %[count edges mid loc] = histcn(VABS(:),vgrid);
%             volume_of_shell = count*d3v;
            % sum(fd3v)
            % units: [f][v3] = l-3 v-3 v3 = l-3 (# per volume)
            % So this is number of particles per volume = density
            [accum_f edges mid loc] = histcn(VABS(:),v_abs_grid,'AccumData',f.f(:)*d3v);
            % sum(fvd3v)
            % units [f][v][v3] = l-3 v-3 v v3 = vl-3 = l s-1 l-1 = s-1 l-2 (# per area and time)
            [accum_fv edges mid loc] = histcn(VABS(:),v_abs_grid,'AccumData',f.f(:).*VABS(:)*d3v);            
            % sum(fvEd3v)
            % units [f][v][E][v3] = l-3 v-3 v eV v3 = v eV l-3 = l s-1 eV l-1 = eV s-1 l-2 (energy per area and time)
            [accum_fvE edges mid loc] = histcn(VABS(:),v_abs_grid,'AccumData',f.f(:).*VABS(:).*ENERGY(:)*d3v);
            % Since we used equally spaced v's, the new phase space volumes
            % (which are shells), become increasingly big further out
            % (larger vabs): d3v = v2*dv*dOmega. Omega is the solid angle
            % and is the same (4*pi for a sphere) for all new bins since we 
            % only binned by |v|.z.
%             fout.volume = volume_of_shell;            
%             fout.vedges = vgrid;
%             fout.dv = dv;
%             fout.Eedges = Egrid;
%             fout.dE = Edelta;
%             fout.vcenter = vcenter;
%             fout.Ecenter = Ecenter;
%             fout.accum_fd3v = accum_f;
            % units 
            % [solid angle] = sr (steradians)
            % [dv] = l s-1
            % [E] = eV
%             fout.dfdv = accum_f/(4*pi)./dv; % m-3/(m/s)= s/m4
%             fout.dfdE = accum_f/(4*pi)./Edelta; % s/(m^4 sr eV)
%             fout.dfvdv = accum_fv/(4*pi)./dv; % (s-1 m-2)/(m/s) = 1/m3
%             fout.dfvdE = accum_fv/(4*pi)./Edelta; % 1/(sr eV) = 1/(s m2 sr eV)
%             fout.dfvEdv = accum_fvE/(4*pi)./dv; % 1/(m/s)
%             fout.dfvEdE = accum_fvE/(4*pi)./Edelta; % (eV s-1 m-2)/(sr eV) = eV/(s m2 sr eV)
            f_vabs_arr(id_count,:) = accum_f/(4*pi)./dv_abs; % m-3/(m/s)= s/m4
            dpf_E_arr(id_count,:) = accum_fv/(4*pi)./Edelta'; % 1/(sr eV) = 1/(s m2 sr eV)
            def_E_arr(id_count,:) = accum_fvE/(4*pi)./Edelta'; % (eV s-1 m-2)/(sr eV) = eV/(s m2 sr eV)
          end
          if any([doPar doPitch])
            if doInterp
              vv = linspace(f.v(1),f.v(end),nInterp);
              [VXm,VYm,VZm] = meshgrid(f.v,f.v,f.v);
              [VXi,VYi,VZi] = meshgrid(vv,vv,vv);
              ff = interp3(VXm,VYm,VZm,f.f,VXi,VYi,VZi);                            
            else
              vv = f.v;
              ff = f.f;
            end
            
            dv = vv(2)-vv(1);                 
            %dv = f.v(2)-f.v(1);
            d3v = dv.^3;
            
            Bx = B{1}; Bx = Bx(id_count);
            By = B{2}; By = By(id_count);
            Bz = B{3}; Bz = Bz(id_count);
            [VX,VY,VZ] = ndgrid(vv,vv,vv);
            Babs = sqrt(Bx.^2 + By.^2 + Bz.^2);
            bx = Bx./Babs;
            by = By./Babs;
            bz = Bz./Babs;
            
            
            
%             
            
            VPAR = VX(:).*bx + VY(:).*by + VZ(:).*bz;
            VABS = sqrt(VX(:).^2 + VY(:).^2 + VZ(:).^2);
            E = VABS.^2/2; % include mass later
            PITCH = acosd(VPAR./VABS);
            
            v_par_grid = linspace(min(VPAR(:)),max(VPAR(:)),n_vpar+1);
            v_par_grid = [vv(1)-0.5*dv_old vv+0.5*dv_old];
            v_par_grid = linspace(vv(1),vv(end),n_vpar+1);
            dv_par = v_par_grid(2)-v_par_grid(1);
            v_par_center = v_par_grid(1:(end-1)) + 0.5*dv_par;
            [accum_fvpar edges mid loc] = histcn(VPAR(:),v_par_grid,'AccumData',ff(:)*d3v);
            f_vpar_arr(id_count,:) = accum_fvpar/dv_par; % m-3/(m/s)= s/m4
            
            % Pitch angle distributions
            % Pitch angle grid
            pitch_grid = linspace(0,180,n_pitch+1);
            d_pitch = pitch_grid(2)-pitch_grid(1);
            pitch_center = pitch_grid(1:(end-1)) + d_pitch;
            % Velocity grid
            v_pitch_grid = linspace(0,1.001*max(VABS(:)),n_v_pitch+1);
            dv_pitch = v_pitch_grid(2)-v_pitch_grid(1);
            v_pitch_center = v_pitch_grid(1:(end-1)) + 0.5*dv_pitch;            
            % Distribution of (vabs,pitch)
            [accum_fpitch edges mid loc] = histcn([VABS(:), PITCH(:)],v_pitch_grid,pitch_grid,'AccumData',ff(:)*d3v);            
            
            % Energy grid            
            E_grid = logspace(log10(0.999*min(E(:))),log10(0.001*max(E(:))),n_E_pitch+1);
            dE = diff(E_grid);
            E_center = E_grid(1:(end-1)) + dE;
            % Distribution of (E,pitch)
            [accum_fEpitch edges mid loc] = histcn([E(:), PITCH(:)],E_grid,pitch_grid,'AccumData',ff(:)*d3v); % 1/vol 
            
            % need to divide this by the opening solid angle/volume
            solang = -(cosd(pitch_grid(2:end)) - cosd(pitch_grid(1:end-1)))*2*pi;
            % velocity
            dvel = v_pitch_grid(2:end).^3/3 - v_pitch_grid(1:end-1).^3/3;
            [SA,DVEL] = meshgrid(solang,dvel);
            DVOL = SA.*DVEL;
            % energy
            dvel = sqrt(2*E_grid(2:end)).^3/3 - sqrt(2*E_grid(1:end-1)).^3/3;
            [SA,DE] = meshgrid(solang,dE);
            DVOLE = SA.*DE;
            
            f_pitch_arr(id_count,:,:) = accum_fpitch./DVOL; % m-3/(sa m / s)= s/m4
            f_pitchE_arr(id_count,:,:) = accum_fEpitch./DVOLE; % m-3/(sr*E)= s/m4
          end
        end 
        
        %[x_arr_sorted,i_sorted] = sort(x_arr);
        i_sorted = 1:numel(x_arr);
        [~,i_sorted] = sort(x_arr);
        
        fout(itime).t = t_arr(i_sorted);          
        fout(itime).x = x_arr(i_sorted);
        fout(itime).z = z_arr(i_sorted);
        fout(itime).v = vaxes;
        fout(itime).fvx = f_vx_arr(i_sorted,:);
        fout(itime).fvy = f_vy_arr(i_sorted,:);
        fout(itime).fvz = f_vz_arr(i_sorted,:);
        if any([doVabs doE]) 
          fout(itime).vabs_center = v_abs_center;
          fout(itime).vabs_edges = v_abs_grid;
          fout(itime).fvabs = f_vabs_arr(i_sorted,:);
          fout(itime).fdpfE = dpf_E_arr(i_sorted,:);
          fout(itime).fdefE = def_E_arr(i_sorted,:);
        end
        if any([doPar doPitch])
          fout(itime).vpar_center = v_par_center;
          fout(itime).vpar_edges = v_par_grid;
          fout(itime).fvpar = f_vpar_arr(i_sorted,:);
          fout(itime).vpitch = v_pitch_grid;
          fout(itime).pitch_edges = pitch_grid;
          fout(itime).pitch_center = pitch_center;
          fout(itime).fpitch = f_pitch_arr(i_sorted,:,:);
          fout(itime).Epitch_edges = E_grid;
          fout(itime).Epitch_center = E_center;
          fout(itime).fpitchE = f_pitchE_arr(i_sorted,:,:);
        end
        if 0
        fout(itime).vabs_edges = vaxes_pos;
        [X,V] = meshgrid(fout(t_count,id).x,vaxes_pos);
        fout(itime).vabs_edges = vaxes_pos;
        fout(itime).vabs_mat = V; % multiply f_vabs_mean_arr(i_sorted,:) with vabs_mat.^4/2 to get DEF        
        fout(itime).fvabssum = f_vabs_sum_arr(i_sorted,:);
        fout(itime).def = f_vabs_mean_arr(i_sorted,:).*V'.^4/2;
        fout(itime).fvabsmean = f_vabs_mean_arr(i_sorted,:);
        end
        %toc        
      end
       
      varargout{1} = fout;
    end
    % Density
    % Flux
    % Velocity
    % Current
    % Stress tensor
    % Pressure    
    % Sets of moments
    
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
  end
  
  methods % set and get
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
    function obj = set.xe1(obj,value)
      obj.xe1_ = value;
    end
    function obj = set.ze1(obj,value)
      obj.ze1_ = value;
    end
    function obj = set.xi1(obj,value)
      obj.xi1_ = value;
    end
    function obj = set.zi1(obj,value)
      obj.zi1_ = value;
    end
    function obj = set.xe2(obj,value)
      obj.xe2_ = value;
    end
    function obj = set.ze2(obj,value)
      obj.ze2_ = value;
    end
    function obj = set.xi2(obj,value)
      obj.xi2_ = value;
    end
    function obj = set.zi2(obj,value)
      obj.zi2_ = value;
    end
    function obj = set.dxi(obj,value)
      obj.dxi_ = value;
    end
    function obj = set.dzi(obj,value)
      obj.dzi_ = value;
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
    function obj = set.id(obj,value)
      obj.id_ = value;
    end
    function obj = set.dists(obj,value)
      obj.dists_ = value;
    end
    function obj = set.tags(obj,value)
      obj.tags_ = value;
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
    function value = get.xe1(obj)
      value = obj.xe1_;
    end
    function value = get.ze1(obj)
      value = obj.ze1_;
    end
    function value = get.xi1(obj)
      value = obj.xi1_;
    end
    function value = get.zi1(obj)
      value = obj.zi1_;
    end
    function value = get.xe2(obj)
      value = obj.xe2_;
    end
    function value = get.ze2(obj)
      value = obj.ze2_;
    end
    function value = get.xi2(obj)
      value = obj.xi2_;
    end
    function value = get.zi2(obj)
      value = obj.zi2_;
    end
    function value = get.dxi(obj)
      value = obj.dxi_;
    end
    function value = get.dzi(obj)
      value = obj.dzi_;
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
    function value = get.id(obj)
      value = obj.id_;
    end
    function value = get.dists(obj)
      value = obj.dists_;
    end
    function value = get.tags(obj)
      value = obj.tags_;
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
    function out = ind_from_val(var,value,varargin)
      % method is the same for xlim, zlim ilim, i, twpelim, twcilim
      
      % Defaults
      doExact = 1;
      doClosest = 0;
      nClosest = 1; % only the closest index
           
      
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
          otherwise
            warning(sprintf('Input ''%s'' not recognized.',args{1}))
            args = args(l+1:end);
        end        
        if isempty(args), break, end    
      end
      
      
      % Find indices
      if doClosest
        inds = find()
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
    function out = get_dists(obj,species)
      % get iterations
      iterations = obj.iteration;
      nIter = obj.length;
      % initialize matrix
      data = nan([obj.get_gridsize,nIter]);
      for iIter = 1:nIter
        iter = iterations(iIter);
        str_iter = sprintf('%010.0f',iter);
        data_tmp = h5read(obj.file,...
           ['/data/' str_iter '/' field],...
           [obj.grid{1}(1) obj.grid{2}(1)],... % start indices
           [numel(obj.grid{1}) numel(obj.grid{2})]); % number of counts
        %disp(sprintf('Reading %s: [%g %g] datapoints starting at [%g %g]',field,numel(obj.grid{1}),numel(obj.grid{2}),obj.grid{1}(1),obj.grid{2}(1)))
        data(:,:,iIter) = data_tmp;
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