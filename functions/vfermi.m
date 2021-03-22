function vafter = vfermi(vbefore,vdf)
% out = vfermi(vbefore,vdf)

for iv = 1:numel(vbefore)
  if vbefore(iv) >= vdf
    vafter(iv) = vbefore(iv);
  else
    vafter(iv) = -vbefore(iv) + 2*vdf;
  end
end