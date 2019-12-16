classdef ObjectArray
   properties
      Value
   end
   methods
      function obj = ObjectArray(F)
         if nargin ~= 0
            m = size(F,1);
            n = size(F,2);
            obj(m,n) = obj;
            for i = 1:m
               for j = 1:n
                  obj(i,j).Value = F(i,j);
               end
            end
         end
      end
   end
end