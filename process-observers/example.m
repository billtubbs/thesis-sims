
function [] = example(varargin)
    for n = 1:2:length(varargin)    
        disp(['Argument ' varargin{n} ' = ' num2str(varargin{n+1})])
    end
