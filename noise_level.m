% Check noise level: std/mean, within given box
varstrs = {'pe1.scalar','pe2.scalar','pi1.scalar','pi2.scalar','E.abs'};
nvars = numel(varstrs);

val_mean = nan(nvars,1);
val_std = nan(nvars,1);

xlim = [20 25];
zlim = [20 25];

ipx1 = find(x>xlim(1),1,'first');
ipx2 = find(x<xlim(2),1,'last');
ipz1 = find(z>zlim(1),1,'first');
ipz2 = find(z<zlim(2),1,'last');
ipx = ipx1:1:ipx2;
ipz = ipz1:1:ipz2;


for ivar = 1:nvars
  varstr = varstrs{ivar};
  variable = eval(varstr);  
  variable_box = variable(ipx,ipz);
  tmp_mean = mean(variable_box(:));
  tmp_std = std(variable_box(:));
  
  val_mean(ivar) = tmp_mean;
  val_std(ivar) = tmp_std;
end