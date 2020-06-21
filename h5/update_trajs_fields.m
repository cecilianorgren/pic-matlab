%tr04 = PICTraj('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5');
tr = tr04.pass('mass',[0.5 1.5]); % dont care about electrons for now
%tr = tr04;
%tr = tr04;
%tr = tr04(find([tr04.id]==723));
isstillnan = struct([]);
for itr = 627:tr.ntr%1%727%627:tr.ntr % 166, 331
  %disp(sprintf('tr = %g/%g',itr,tr.ntr))
  % It will probably be faster to just update the nan indices
  tic
  inanEx = find(isnan(tr(itr).Ex));
  inanEy = find(isnan(tr(itr).Ey));
  inanEz = find(isnan(tr(itr).Ez));
  inanBx = find(isnan(tr(itr).Bx));
  inanBy = find(isnan(tr(itr).By));
  inanBz = find(isnan(tr(itr).Bz));
  inan = unique([inanEx;inanEy;inanEz;inanBx;inanBy;inanBz]);
  %inan = 1:tr(itr).length; % redo all, because i had the wrong z before; up until 166
  disp(sprintf('tr = %g/%g, nnan = %g',itr,tr.ntr,numel(inan)))
  
  if not(isempty(inan))
  %for in = 1:numel(inan)
    [repEx,repEy,repEz,repBx,repBy,repBz] = df04.interp_EB3(tr(itr).x(inan),tr(itr).z(inan),tr(itr).t(inan));
  %end
    Ex = tr(itr).Ex; Ex(inan) = repEx;
    Ey = tr(itr).Ey; Ey(inan) = repEy;
    Ez = tr(itr).Ez; Ez(inan) = repEz;
    Bx = tr(itr).Bx; Bx(inan) = repBx;
    By = tr(itr).By; By(inan) = repBy;
    Bz = tr(itr).Bz; Bz(inan) = repBz;
    plot(tr(itr).t,Ez,tr(itr).t,tr(itr).Ez)
    drawnow;
    inanEx = find(isnan(Ex));
    inanEy = find(isnan(Ey));
    inanEz = find(isnan(Ez));
    inanBx = find(isnan(Bx));
    inanBy = find(isnan(By));
    inanBz = find(isnan(Bz));
    inan = unique([inanEx;inanEy;inanEz;inanBx;inanBy;inanBz]);
    if not(isempty(inan))
      disp(sprintf('id = %g, nnan = %g',tr(itr).id,numel(inan)))
      isstillnan(end+1).id = tr(itr).id;
      isstillnan(end+1).nnan = numel(inan);
    end

    if 1
      h5write_trajs_update('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5',tr(itr).id,'Ex',Ex)
      h5write_trajs_update('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5',tr(itr).id,'Ey',Ey)
      h5write_trajs_update('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5',tr(itr).id,'Ez',Ez)
      h5write_trajs_update('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5',tr(itr).id,'Bx',Bx)
      h5write_trajs_update('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5',tr(itr).id,'By',By)
      h5write_trajs_update('/Volumes/Fountain/Data/PIC/df_cold_protons_n04/data_h5/trajectories.h5',tr(itr).id,'Bz',Bz)
    end
  end
  toc
end
disp('Done.')