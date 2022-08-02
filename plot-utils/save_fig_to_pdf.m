function save_fig_to_pdf(filename, h)
% save_fig_to_pdf(filename, h)
% Saves a pdf document of the current figure with the paper
% size adjusted to the size of the figure.
    if nargin < 2
        h = gcf;
    end
    % This existing function by Gabe Hoffmann works better
    save2pdf(filename,h)
end