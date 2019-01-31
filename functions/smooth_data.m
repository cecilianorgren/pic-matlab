function out = smooth_data(data,varargin)

new_data = smooth(data,varargin{:});
new_data = reshape(new_data,size(data));

out = new_data;