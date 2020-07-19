
str = 'log10(Bx./Bz)';
str = 'Bx./Bz';
met = methods(pic);
found = cellfun(@(s) strfind(str,s),met,'UniformOutput',false);
%found_ = cellfun(@(s) find(strfind(str,s)),met);
found_ = cellfun(@(s) contains(str,s),met);
list_fields = met(cellfun(@(s) contains(str,s),met))
%operators = split(str,met);
%operators(find(cellfun(@(s) isempty(s),operators))) = []; % remove empty

%strfind(str,spl)