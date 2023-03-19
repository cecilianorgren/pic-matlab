function colors = pic_colors(colorOption)
% MMS_COLORS 
%   colors = mms_colors('1234'); % mms1, mms2, mms3, mms4
%   colors = mms_colors('xyza'); % x, y, z, abs
%   colors = mms_colors('12yzab'); % mms1, mms2, y, z, abs, unknown (=yellow)
path_pic = fileparts(which('pic_startup.m'));

switch colorOption
  case 'matlab'
    colors = [0.0000    0.4470    0.7410;...
              0.8500    0.3250    0.0980;...
              0.9290    0.6940    0.1250;...
              0.4940    0.1840    0.5560;...
              0.4660    0.6740    0.1880;...
              0.3010    0.7450    0.9330;...
              0.6350    0.0780    0.1840];
  case 'extended'
    colors = [0.00  0.00  1.00;...
              0.00  0.50  0.00;...
              1.00  0.00  0.00;... 
              0.00  0.75  0.75;...
              0.75  0.00  0.75;...
              0.75  0.75  0.00;...
              0.25  0.25  0.25;...
              0.75  0.25  0.25;...
              0.95  0.95  0.00;... 
              0.25  0.25  0.75;...
              0.75  0.75  0.75;...
              0.00  1.00  0.00;... 
              0.76  0.57  0.17;...
              0.54  0.63  0.22;...
              0.34  0.57  0.92;...
              1.00  0.10  0.60;...
              0.88  0.75  0.73;...
              0.10  0.49  0.47;...
              0.66  0.34  0.65;...
              0.99  0.41  0.23];
            
  case {'bluegrayred'}
      rr = interp1([1 64 128 192 256],[0.0  0.5 0.9 1.0 0.75],1:256);
      gg = interp1([1 64 128 192 256],[0.0  0.5 0.9 0.5 0.00],1:256);
      bb = interp1([1 64 128 192 256],[0.75 1.0 0.9 0.5 0.00],1:256);
      colors = [rr' gg' bb'];      
  case 'blue_red'
    out = load([path_pic filesep 'colormaps' filesep 'blue_red.mat']);
    colors = out.cmap;
  case 'blue_white'
    out = load([path_pic filesep 'colormaps' filesep 'blue_white.mat']);
    colors = out.cmap;
  case 'islands'
    out = load([path_pic filesep 'colormaps' filesep 'white_jet.mat']);
    colors = out.cmap;    
  case 'blue_gray_red' % fancy-schmancy
    c = [0,0,1;...
        0.92,0.92,0.92;...
        1,0,0];
    colors = interp1(linspace(1,64,size(c,1)),c,1:64);  
  case 'waterfall' % fancy-schmancy
    c = [55,137,187;...
        106,193,165;...
        172,220,166;...
        230,244,157;...
        255,254,194;...
        253,223,144;...
        251,173,104;...
        242,109,074;...
        211,064,082]/255;
    colors = interp1(linspace(1,64,size(c,1)),c,1:64);  
  case 'waterfall2' % fancy-schmancy
    c = [55,137,187;...
        106,193,165;...
        172,220,166;...
        230,244,157;...
        255,254,194;...
        253,223,144;...
        251,173,104;...
        242,109,074;...
        211,064,082;...
        100,107,137;...
        ...%55,137,187;...
        10 100 150]/255;
    colors = interp1(linspace(1,64,size(c,1)),c,1:64);  
  case 'candy' % fancy-schmancy
    c = [255 255 255;     
        106,193,165;...
        172,220,166;... 
        %255,254,194;...
        230,244,157;...               
        253,223,144;...
        251,173,104;...
        242,109,074;...
        211,064,082;...
        %100,107,137;...
        %55,137,187...
        ]/255;
    colors = interp1(linspace(1,64,size(c,1)),c,1:64);  
  case 'candy3' % fancy-schmancy
    c = [255 255 255;...
        200 200 200;...      
        106,193,165;...
        106,193,165;...
        172,220,166;... 
        %255,254,194;...
        230,244,157;...               
        253,223,144;...
        251,173,104;...
        242,109,074;...
        211,064,082;...
        100,107,137;...
        ...%55,137,187;...
        10 100 150]/255;
    colors = interp1(linspace(1,64,size(c,1)),c,1:64);  
  case 'candy2' % fancy-schmancy with white in the middle
    c = [%230,244,157;...               
        %70,160,85;...
        106,193,165;...
        172,220,166;... 
        253,223,144;...
        243,233,144;...
        %255,254,194;...
        %255 255 255;...
        255 255 255;...
        251,173,104;...
        242,109,074;...
        211,064,082;...        
        100,054,082;...
        %100,107,137;...
        %75,107,187 ... % 55,137,187 ...
        ]/255;
    colors = interp1(linspace(1,101,size(c,1)),c,1:101);
    % add 2 extra white rows to make '0' to make it more visible
    [C,IA,IB] = intersect(colors,[1 1 1],'rows');
    %colors = [colors(1:(IA-1),:); C; C; C; colors((IA+1):end,:)];
  case 'candy4' % fancy-schmancy, with dark red/black on the end, for contrast
    c = [255 255 255;     
        106,193,165;...
        172,220,166;... 
        %255,254,194;...
        230,244,157;...               
        253,223,144;...
        251,173,104;...
        242,109,074;...
        211,064,082;...
        100,054,082;...
        %100,107,137;...
        %55,137,187...
        ]/255;
      colors = interp1(linspace(1,101,size(c,1)),c,1:101);
      % add 2 extra white rows to make '0' to make it more visible
      [C,IA,IB] = intersect(colors,[1 1 1],'rows');
      colors = [colors(1:(IA-1),:); C; C; C; colors((IA+1):end,:)];      
  case 'candy5' % fancy-schmancy, with dark red/black on the end, for contrast
    c = [255 255 255;
        253,223,144;...
        251,173,104;...
        242,109,074;...
        211,064,082;...
        100,054,082;...        
        0.3*[100,054,082];...
        %10,10,10;...
        %100,107,137;...
        %55,137,187...
        ]/255;
      colors = interp1(linspace(1,101,size(c,1)),c,1:101);
      % add 2 extra white rows to make '0' to make it more visible
      [C,IA,IB] = intersect(colors,[1 1 1],'rows');
      colors = [colors(1:(IA-1),:); C; C; C; colors((IA+1):end,:)];   
  case 'thermal' % fancy-schmancy, with dark red/black on the end, for contrast
    c = [255 255 255;
        253,223,144;...
        251,173,104;...
        242,109,074;...
        211,064,082;...
        100,054,082;...        
        0.3*[100,054,082];...
        %10,10,10;...
        %100,107,137;...
        %55,137,187...
        ]/255;
      c = flipdim(c,1);
      colors = interp1(linspace(1,101,size(c,1)),c,1:101);
      % add 2 extra white rows to make '0' to make it more visible
      [C,IA,IB] = intersect(colors,[1 1 1],'rows');
      colors = [colors(1:(IA-1),:); C; C; C; colors((IA+1):end,:)]; 
  case 'cycle'
        c = [0.0000    0.4470    0.7410;...
             0.4660    0.6740    0.1880;...
             0.9290    0.8940    0.1250;...
             0.9290    0.8940    0.1250;...
             0.8500    0.3250    0.0980;...
             0.0000    0.4470    0.7410];
    colors = interp1(linspace(1,64,size(c,1)),c,1:64);  
  case 'cycle2'
        c = [0.2157    0.5373    0.7333;...
             0.5    0.9    0.5;...             
             1.0000    0.9961    0.5608;...
             1.0000    0.9961    0.5608;...
             0.8275    0.2510    0.3216;...
             0.2157    0.5373    0.7333];
    colors = interp1(linspace(1,64,size(c,1)),c,1:64);  
  case 'pasteljet'
    colors = customcolormap([0 .25 .5 .75 1], {'#9d0142','#f66e45','#ffffbb','#65c0ae','#5e4f9f'});
  case 'red-yellow-blue'
    colors = customcolormap(linspace(0,1,11), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd','#def4f9','#abd9e9','#73add2','#4873b5','#313691'});
  case 'red-yellow-green'
    colors = customcolormap(linspace(0,1,11), {'#a60126','#d7302a','#f36e43','#faac5d','#fedf8d','#fcffbf','#d7f08b','#a5d96b','#68bd60','#1a984e','#006936'});
  case 'red-white-blue'
    colors = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});
  case 'orange-white-purple'
    colors = customcolormap(linspace(0,1,11), {'#7f3c0a','#b35807','#e28212','#f9b967','#ffe0b2','#f7f7f5','#d7d9ee','#b3abd2','#8073a9','#562689','#2f004d'});
  case 'purple-white-green'
    colors = customcolormap(linspace(0,1,11), {'#410149','#762a84','#9b6fac','#c1a5cd','#e7d4e8','#faf6f7','#d7f1d6','#a6db9d','#5aae60','#1c7735','#014419'});
  case 'pink-white-green'
    colors = customcolormap(linspace(0,1,11), {'#860454','#c51b7c','#dc75ab','#f0b7da','#ffdeef','#f8f7f7','#e5f4d9','#b9e084','#7fbc42','#4d921e','#276418'});
  case 'brown-white-pool'
    colors = customcolormap(linspace(0,1,11), {'#523107','#523107','#bf812c','#e2c17e','#f3e9c4','#f6f4f4','#cae9e3','#81cdc1','#379692','#01665e','#003d2e'});
  otherwise % abcd,xyz,1234,g
    colors = [];
    for ii = 1:numel(colorOption)
      switch colorOption(ii)
        case 'a' % hot ions
          newColor = [0.8500    0.3250    0.0980];
          %newColor = [0.4940    0.1840    0.5560];
          %newColor = [147    112 219]/255;  
        case 'b' % hot electrons
          newColor = [0.9290    0.6940    0.1250];
          newColor = [219    112    147]/255;
          newColor = [255    182    193]/255;
        case 'c' % cold ions
          newColor = [0.0000    0.4470    0.7410];
          %newColor = [147    112 219]/255;          
        case 'd' % cold electrons
          newColor = [0.3010    0.7450    0.9330];
          %newColor = [ 0.4660    0.6740    0.1880];          
        case '1' % mms1, black
          newColor = [0 0 0];
        case '2' % mms2, red
          newColor = [1 0.2 0];
        case '3' % mms3, bluish green
          newColor = [0 0.8 0.0];
        case '4' % mms4, sky blue
          newColor = [0.1 0.4 1];
        case 'x' % x
          newColor = [0.1 0.4 1];
        case 'y' % y
          newColor = [0 0.7 0.0];
        case 'z' % z
          newColor = [1 0.2 0];
        case 'g' % absolute value
          newColor = [0.7 0.7 0.7]; 
        case '0'
           newColor = [   0.9290    0.6940    0.1250];
        case 'm' % purple
          newColor = [ 0.5540    0.2340    0.6060];
        case 'p' % light purple
          newColor = [216 191 216]/255;  
        otherwise 
          irf.log('warning',['Can''t recognize input ''' colorOption(ii) '''. Inserting a nice yellow instead.']);
          newColor = [0.95 0.7 0.0]; 
      end  
      colors(ii,1:3) = newColor;  
    end
end