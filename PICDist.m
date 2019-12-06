 classdef PICDist
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
    id_
    dists_
    
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
    id    
    dists
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
            
      %obj.charge = obj.get_charge;
      %obj.mass = obj.get_mass;
      %uniqueMass = sort(unique(obj.mass));
      %obj.mime = uniqueMass(2)/uniqueMass(1); % second lightest/lightest
      %obj.teti = h5read(h5filePath,'/simulation_information/teti');
      %obj.wpewce = h5read(h5filePath,'/simulation_information/wpewce');
      
      obj.iteration = get_iterations(obj);    
      obj.it = 1:numel(obj.iteration);
      obj.dists_ = obj.get_distlist;
      [xi,zi] = obj.get_locs;
      obj.xi = xi;
      obj.zi = zi;
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
          obj.xi_ = builtin('subsref',obj.xi_,s);
          obj.zi_ = builtin('subsref',obj.zi_,s);
          obj.it_ = builtin('subsref',obj.it_,s);
         
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
    function obj = update_inds(obj,inds)
      % xi
      % zi
      % dists
      % indices
      for it = 1:obj.nt        
        obj.xi_{it} = obj.xi_{it}(inds{it});
        obj.zi_{it} = obj.zi_{it}(inds{it});
        obj.dists_{it} = obj.dists_{it}(inds{it});
        obj.indices_{it} = obj.indices_{it}(inds{it});
      end
      
      
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
    function varargout = plot_map(obj,iSpecies,sumdim)
            
      fontsize = 7;
      %doBDir = 1;
      ticks = -15:1:15;
      
      xlim = [min([obj.xi_{1}{:}]) max([obj.xi_{1}{:}])];
      zlim = [min([obj.zi_{1}{:}]) max([obj.zi_{1}{:}])];      
      fig_position = get(0,'screensize'); %[1 330 2555 1015]; 
      fig = figure;
      fig.Position = fig_position; 

      idist = 0;
      %tic
      h = gobjects(0);
      hleg = gobjects(0);
      for id = obj.indices{1}
        idist = idist + 1;      

        % Load data          
        f = obj.fxyz(1,idist,iSpecies,sumdim); % sum over 3rd dim
        %if not(all(xlo>xlim(1) && xhi<xlim(2) && zlo>zlim(1) && zhi<zlim(2)))
        %  disp([sprintf('%.2f %.2f %.2f %.2f outside of box',xlo,xhi,zlo,zhi)])
        %  continue
        %end
        xloc = (f.x(1)+f.x(2))/2;
        zloc = (f.z(1)+f.z(2))/2;
        %disp(sprintf('%.2f ',f.x(1),f.x(2),f.z(1),f.z(2)))
        axes_position = [(f.x(1)-xlim(1))/diff(xlim) ...
                         (f.z(1)-zlim(1))/diff(zlim) ...
                         (f.x(2)-f.x(1))/diff(xlim) ...
                         (f.z(2)-f.z(1))/diff(zlim)];

        pause(0.1)  

        if 1 % xz plane
          nrows = 1;
          ncols = 3;
          npanels = nrows*ncols;

          hca = subplot('Position',axes_position);
          h(end+1) = hca;
          hca = gca;

            if 0 % fxy
              hca = h(isub); isub = isub + 1;
              imagesc(hca,vx(:,ispecies),vy(:,ispecies),squeeze(fxy(:,:,ispecies))')
              hca.XLabel.String = 'vx';
              hca.YLabel.String = 'vy';
              hcb = colorbar('peer',hca);
              hcb.YLabel.String = 'f';
              hca.Title.String = strtitle;
            end
            if 0 % fxy
              imagesc(hca,vx(:,ispecies),vy(:,ispecies),squeeze(fxy(:,:,ispecies))')  
              hca.XLabel.String = '';
              hca.YLabel.String = '';
              %hcb = colorbar('peer',hca);
              %hcb.YLabel.String = 'f';
              %hca.Title.String = strtitle;
              hca.YDir = 'normal';        
              hca.XLim = vlim(ispecies)*[-1 1];
              hca.YLim = vlim(ispecies)*[-1 1];
              hca.XGrid = 'on';
              hca.YGrid = 'on';
              hca.XTick = ticks;
              hca.YTick = ticks;
              %if 
              hca.XTickLabel = [];
              hca.YTickLabel = [];
              hca.Box = 'on';
              colormap(hca,pic_colors('candy'))
              irf_legend(hca,{sprintf('x=%.1f, z=%.1f',xloc,zloc);sprintf('B=[%.2f,%.2f,%.2f]',Bloc.x,Bloc.y,Bloc.z)},[0.01 0.99],'color',[0 0 0],'fontsize',fontsize)

              if doBDir           
                line_slope = (Bloc.y/Bloc.x);
                xx = min(hca.XLim(2)*[1 1/abs(line_slope)])*[-1 1];
                hold(hca,'on')                
                hBline = plot(hca,xx,xx*line_slope,'linewidth',0.5,'color',[0.5 0.5 0.5]);
                hold(hca,'off')
              end
            end
            if 1 % fxz
              imagesc(hca,f.v,f.v,f.f')
              hca.XLabel.String = '';
              hca.YLabel.String = '';
              hca.YDir = 'normal';        
              %hca.XLim = vlim(ispecies)*[-1 1];
              %hca.YLim = vlim(ispecies)*[-1 1];
              hca.XGrid = 'on';
              hca.YGrid = 'on';
              hca.XTick = ticks;
              hca.YTick = ticks;              
              hca.XTickLabel = [];
              hca.YTickLabel = [];
              hca.Box = 'on';
              colormap(hca,pic_colors('candy'))
              %irf_legend(hca,{sprintf('x=%.1f, z=%.1f',xloc,zloc);sprintf('B=[%.2f,%.2f,%.2f]',Bloc.x,Bloc.y,Bloc.z)},[0.01 0.99],'color',[0 0 0],'fontsize',9)
              hleg_ = irf_legend(hca,{sprintf('x=%.1f, z=%.1f',xloc,zloc)},[0.01 0.99],'color',[0 0 0],'fontsize',fontsize);
              hleg(end+1) = hleg_;
              if 0%doBDir           
                line_slope = (Bloc.z/Bloc.x);
                xx = min(hca.XLim(2)*[1 1/abs(line_slope)])*[-1 1];
                hold(hca,'on')                
                hBline = plot(hca,xx,xx*line_slope,'linewidth',0.5,'color',[0.5 0.5 0.5]);
                hold(hca,'off')
              end
            end
            if 0 % fzy
              hca = h(isub); isub = isub + 1;
              imagesc(hca,vy(:,ispecies),vz(:,ispecies),squeeze(fyz(:,:,ispecies))')  
              hca.XLabel.String = 'vz';
              hca.YLabel.String = 'vy';
              hcb = colorbar('peer',hca);
              hcb.YLabel.String = 'f';
              hca.Title.String = strtitle;
            end
            %print('-dpng','-r200',[savedir_root sub_dir '/' strprint '.png']);
            drawnow
            %pause(1)
        end
      end
        %toc      
      
      varargout{1}.h = h;
      varargout{1}.hleg = hleg;
      
    end
              
    % Data analysis routines, time derivatives, interpolation, etc.
    function out = get_peaks(obj,nPeaks,spacingPeaks,iSpecies)
      nTimes = obj.nt;
      nDists = obj.nd;
      fpeaks = struct;
      for it = 1:nTimes
        for id = 1:nDists{it}
          f = obj.f(it,id,iSpecies);
          ftmp = f.f;
          for iPeak = 1:nPeaks
            [val,ind] = max(ftmp(:));
            [ix,iy,iz] = ind2sub(size(ftmp),ind);
            ix_ = ix+[-spacingPeaks:spacingPeaks];
            iy_ = iy+[-spacingPeaks:spacingPeaks];
            iz_ = iz+[-spacingPeaks:spacingPeaks];
            ftmp(ix_,iy_,iz_) = NaN;
            fpeaks(iPeak,id,it).vx = f.v(ix);
            fpeaks(iPeak,id,it).vy = f.v(iy);
            fpeaks(iPeak,id,it).vz = f.v(iz);
            fpeaks(iPeak,id,it).x = mean(f.x);
            fpeaks(iPeak,id,it).y = 0;
            fpeaks(iPeak,id,it).z = mean(f.z);
            fpeaks(iPeak,id,it).f = val;
            
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
    
    function [x,z] = get_locs(obj)
      
      iterations = obj.iteration;
      x = cell(obj.nt,1);
      z = cell(obj.nt,1);
      for iIter = 1:obj.nt
        iter = iterations(iIter);
        str_iter = sprintf('%010.0f',iter);        
        for iDist = 1:numel(obj.dists{iIter})
          dataset_name = ['/data/' str_iter '/' num2str(iDist,'%05.0f') '/fxyz'];
          %locs{iIter}{iDist} = [h5readatt(obj.file,dataset_name,'x') h5readatt(obj.file,dataset_name,'z')];
          x{iIter}{iDist} = h5readatt(obj.file,dataset_name,'x')';
          z{iIter}{iDist} = h5readatt(obj.file,dataset_name,'z')';
        end
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
    
    % dEF - differential energy flux
    function varargout = dEF(obj,it,id,iss)
      % out = dEF(obj,it,id,iss)
      f = obj.fxyz(it,id,iss);
      nEnergy = 50;
      %energy = logspace(-3,log10(max(f.v.^2)),nEnergy);
      [VX,VY,VZ] = ndgrid(f.v,f.v,f.v);
      VV2 = VX.^2 + VY.^2 + VZ.^2;
      fvv = f.f.*VV2;
      ENERGY = VV2/2;
      energy_edges = logspace(-3,log10(1.0*max(f.v.^2/2)),30);
      energy_edges = logspace(-2,log10(0.2*max(f.v.^2/2)),50);
      %ienergy = hist(energy(:,ispecies))
      [N,EDGES,BIN] = histcounts(ENERGY(:),energy_edges);
      f_energy_edges = tocolumn(EDGES);
      f_energy_centers = tocolumn((EDGES(2:end)+EDGES(1:end-1))*0.5);
      nbins = (numel(energy_edges)-1);
      for ibin = 1:nbins
        ind_bin = find(BIN==ibin);      
        f_dist_tmp = f.f(ind_bin);
        f_dist_mean(ibin,1) = mean(f_dist_tmp);
        f_dist_sum(ibin,1) = sum(f_dist_tmp);
      end
      if nargout == 1
        varargout{1}.energy = f_energy_centers;
        varargout{1}.energy_edges = f_energy_edges;
        varargout{1}.dEF = f_dist_sum;
        varargout{1}.x = f.x;
        varargout{1}.z = f.z;
      elseif nargout == 2        
        varargout{2}.energy_edges = f_energy_edges;
        varargout{3}.dEF = f_dist_sum;
      elseif nargout == 3
        varargout{1}.energy = f_energy_centers;
        varargout{2}.energy_edges = f_energy_edges;
        varargout{3}.dEF = f_dist_sum;
      end
        
%       f_dist_all.distnumber(idist,1) = distnumber;
%       f_dist_all.x(idist,1) = xc;
%       f_dist_all.z(idist,1) = zc;
%       f_dist_all.f_energy_centers(idist,ispecies,:) = f_energy_centers{ispecies};
%       f_dist_all.f_energy_edges(idist,ispecies,:) = f_energy_edges{ispecies};
%       f_dist_all.f_dist_sum(idist,ispecies,:) = f_dist_sum{ispecies};
%       f_dist_all.f_dist_mean(idist,ispecies,:) = f_dist_mean{ispecies};
    
      
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
    function obj = set.id(obj,value)
      obj.id_ = value;
    end
    function obj = set.dists(obj,value)
      obj.dists_ = value;
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
    function value = get.id(obj)
      value = obj.id_;
    end
    function value = get.dists(obj)
      value = obj.dists_;
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