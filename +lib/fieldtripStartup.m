function fieldtripStartup(ftPath,varargin)

restoredefaultpath
if ((nargin>0) && ~isempty(ftPath))
    addpath(ftPath);
    ft_defaults
end

if (nargin<2)
    return;
else
    for i = 2:nargin
        addpath(varargin{i-1});
    end
end

clc;

end