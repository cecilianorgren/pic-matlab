%% EKSPORTER TIL PNG
function lagenpng_keepsize(namn)
% figuresize(42.0/2,29.7/2.7,'centimeters')
% figuresize(42.0/1.2,29.7/2.7,'centimeters')

set(gcf, 'Color', 'w');
str1=[num2str(namn),'.png'];

export_fig(str1,'-png','-q101','-r250')

