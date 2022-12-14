mime = 100;

xmax_de = 512/2;
xmin_de = 0;
zmin_de = -128/2;
zmax_de = 128/2;

xmin_di = xmin_de/sqrt(mime);
xmax_di = xmax_de/sqrt(mime);
zmin_di = xmin_de/sqrt(mime);
zmax_di = xmax_de/sqrt(mime);

dx = 0.25;
dz = 0.25;

z_center = (zmax_di-zmin_di)/2;
x_center = (xmax_di-xmin_di)/2;

x_box = x_center + [-0.5 0 0.5];%(-3:(2*dx):3);
z_box = z_center + (-3:(2*dz):3);

[X,Z] = ndgrid(x_box,z_box);

all_boxes = [X(:)-dx X(:)+dx Z(:)-dz Z(:)+dz];



