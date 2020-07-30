pic = nobg;
twpe = [9000];
xlim = 103+0.5*[-1 1];
zlim = [0 15];
pic = pic.xlim(xlim).zlim(zlim).twpelim(twpe,'exact');

str = 'log10(Bx./Bz)';
str = 'n(1)-Ey+tzz([3 5])';
eval_str = str;
met = methods(pic);
found = cellfun(@(s) strfind(str,s),met,'UniformOutput',false);
%found_ = cellfun(@(s) find(strfind(str,s)),met);
found_ = cellfun(@(s) contains(str,s),met);
list_fields = met(cellfun(@(s) contains(str,s),met));
list_fields = regexpi(str,{'[a-z_A-Z]*'},'match');
list_fields = intersect(list_fields{1},met);

varcount = 0;
for ifield = 1:numel(list_fields)
  nchar = numel(list_fields{ifield});
  str_loc = strfind(str,list_fields{ifield});
  for istr_ = 1:numel(str_loc)
    istr = str_loc(istr_);
    varcount = varcount + 1; 
    if strfind(str(istr+nchar),'(')
      ind1 = strfind(str(istr+nchar:end),'(');
      ind2 = strfind(str(istr+nchar:end),')');
      %varsplit = regexp(varstrs{ivar}, '(?<var>\w+)\W+(?<ind>\d)\W+','names');
      indstr = str(istr + nchar - 1 + (ind1(1)+1:ind2(1)-1));
      indstr
      tmp_str = [list_fields{ifield} '(' indstr ')'];
      
      var = pic.(list_fields{ifield})(eval(indstr));
    else
      var = pic.(list_fields{ifield});
      tmp_str = list_fields{ifield};
    end 
    eval_str = strrep(eval_str,tmp_str,['allvars{' num2str(varcount) '}']);
    allvars{varcount} = var;
  end
end
out = eval(eval_str);

       

%operators = split(str,met);
%operators(find(cellfun(@(s) isempty(s),operators))) = []; % remove empty

%strfind(str,spl)