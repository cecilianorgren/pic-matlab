function varargout = smooth_(varargin)

out = smooth(varargin{1},'moving',varargin{2});

varargout{1} = out;