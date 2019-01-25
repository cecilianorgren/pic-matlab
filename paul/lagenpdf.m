%% EKSPORTER TIL PDF
function lagenpdf(namn)

% figuresize(42.0/2,29.7/2,'centimeters')
% figuresize(42.0/3.5,29.7/3.5,'centimeters')
% figuresize(42.0/3,29.7/4.5,'centimeters')

set(gcf, 'Color', 'w');
%export_fig(sprintf(namn),'-pdf','-transparent','-painters')
export_fig(sprintf(namn),'-pdf')