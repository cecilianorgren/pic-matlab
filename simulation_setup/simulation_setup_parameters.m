% Simulation setup
clear sim sim_compact
% Predefined in code

% sim.wpewce = 2;
% sim.dtwpe = 0.5;
% sim.mi = 100;
% sim.me = 1;
% sim.xemin = 0;
% sim.xemax = 2048;
% sim.zemin = -128;
% sim.zemax = 128;
% sim.nx = 6400;
% sim.nz = 1600;
% sim.nsp = 6; % number of species

sim.wpewce = 2;
sim.dtwpe = 0.5;
sim.mi = 100;
sim.me = 1;
sim.xemin = 0;
sim.xemax = 2*2048;
sim.zemin = -300;-3*128;
sim.zemax = 300;3*128;
sim.nx = 2.0*6400;
sim.nz = 2.0*1600;
sim.nsp = 6; % number of species

% Derived
sim.mime = sim.mi/sim.me;
sim.ximin = sim.xemin/sqrt(sim.mime);
sim.ximax = sim.xemax/sqrt(sim.mime);
sim.zimin = sim.zemin/sqrt(sim.mime);
sim.zimax = sim.zemax/sqrt(sim.mime);
sim.dxe = (sim.xemax-sim.xemin)/sim.nx;
sim.dze = (sim.zemax-sim.zemin)/sim.nz;
sim.dxi = (sim.ximax-sim.ximin)/sim.nx;
sim.dzi = (sim.zimax-sim.zimin)/sim.nz;
sim.dtwci = sim.dtwpe/(2*sim.mime);

% Computer setup
% task = proc?
sim.np_per_task = 9e8;%320000000;
sim.n_nodes = 80;
sim.n_tasks_per_node = 1; % every task needs it own copy of the fields
sim.n_cpus_per_task = 32; % these can share memory between them
sim.n_tasks = sim.n_nodes*sim.n_tasks_per_node;
sim.n_cpus = sim.n_tasks*sim.n_cpus_per_task;
sim.tot_num_particles = sim.n_tasks*sim.np_per_task;

sim.n_particle_files = sim.n_tasks;

% Memory
% Each task/proc needs the field files and the particle files
% Both fields and particle quantities are declared as 'real' in globalshf.h.
% 1 byte = 8 bits
% 1 byte = 2^30 GB ?
% 1 byte = 10^6 GB ?
sim.byte_per_GB = 10^9;2^30;

% Grid 6400x3200 gives files of 2.7 GB
% Quantities stored: (6+10*nsp)*nx*nz
% bx,by,bz,ex,ey,ez = 6
% dns,vxs,vys,vzs,pxx,pyy,pzz,pxy,pxz,pyz = 10*nsp
sim.byte_per_entry_fields = 4; % 128 bits
sim.n_entries_field_files = (6+10*sim.nsp)*sim.nx*sim.nz;
sim.size_field_files_bytes = sim.byte_per_entry_fields*sim.n_entries_field_files;
sim.size_field_files_GB = sim.size_field_files_bytes/sim.byte_per_GB;

% np = 320000000 gives files of 6.4 GB
% Quantities stored: np*5
% x,z,vx,vy,vz
sim.byte_per_entry_particles = 4; % 32 bits
sim.n_entries_particle_files = 5*sim.np_per_task; % x, z, vx, vy, vz
sim.size_particle_files_bytes = sim.byte_per_entry_particles*sim.n_entries_particle_files;
sim.size_particle_files_GB = sim.size_particle_files_bytes/sim.byte_per_GB;

% Fram memory (min/max) = 64/512 (512 is for some special nodes)
sim.additional_memory_overhead_percentage = 100;
sim.memory_usage_per_node_GB = sim.n_tasks_per_node*(sim.size_particle_files_GB + sim.size_field_files_GB)*(1 + sim.additional_memory_overhead_percentage/100);
sim.memory_per_node_fram = 64;
sim.memory_per_node_vilje = 32;

% Particles per cell 
% Division of particles between species
sim.n_part(1) = sim.tot_num_particles*(0.5/9);
sim.n_part(2) = sim.tot_num_particles*(0.5/9);
sim.n_part(3) = sim.tot_num_particles*(2/9);
sim.n_part(4) = sim.tot_num_particles*(2/9);
sim.n_part(5) = sim.tot_num_particles*(2/9);
sim.n_part(6) = sim.tot_num_particles*(2/9);
sim.n_part_per_cell_average = round(sim.n_part/sim.nx/sim.nz);
% I'm not sure this meets up with the actual number of particles assigned

% More compact structure with most relevant properties
fields_to_copy = {'mime','wpewce','nx','nz','dtwpe','dtwci',...
  'dxe','dze','dxi','dzi','ximin','ximax','zimin','zimax',...
  'np_per_task','n_nodes',...
  'n_tasks_per_node','n_tasks','n_cpus_per_task','n_cpus',...
  'size_field_files_GB','size_particle_files_GB','additional_memory_overhead_percentage',...
  'memory_usage_per_node_GB',...
  'memory_per_node_fram','n_part_per_cell_average'};
for ifield = 1:numel(fields_to_copy)
  sim_compact.(fields_to_copy{ifield}) = sim.(fields_to_copy{ifield});
end

%sim
sim_compact

%% Start from the other direction, with particles per cell etc.
sim.n_part_per_cell_average = [100 100 400 400 400];

% Simulation setup
clear sim sim_compact
% Predefined in code

% sim.wpewce = 2;
% sim.dtwpe = 0.5;
% sim.mi = 100;
% sim.me = 1;
% sim.xemin = 0;
% sim.xemax = 2048;
% sim.zemin = -128;
% sim.zemax = 128;
% sim.nx = 6400;
% sim.nz = 1600;
% sim.nsp = 6; % number of species

sim.wpewce = 2;
sim.dtwpe = 0.5;
sim.mi = 100;
sim.me = 1;
sim.xemin = 0;
sim.xemax = 2*2048;
sim.zemin = -3*128;
sim.zemax = 3*128;
sim.nx = 2.0*6400;
sim.nz = 2*1600;
sim.nsp = 6; % number of species

% Derived
sim.mime = sim.mi/sim.me;
sim.ximin = sim.xemin/sqrt(sim.mime);
sim.ximax = sim.xemax/sqrt(sim.mime);
sim.zimin = sim.zemin/sqrt(sim.mime);
sim.zimax = sim.zemax/sqrt(sim.mime);
sim.dxe = (sim.xemax-sim.xemin)/sim.nx;
sim.dze = (sim.zemax-sim.zemin)/sim.nz;
sim.dxi = (sim.ximax-sim.ximin)/sim.nx;
sim.dzi = (sim.zimax-sim.zimin)/sim.nz;
sim.dtwci = sim.dtwpe/(2*sim.mime);

% Computer setup
% task = proc?
sim.np_per_task = 9e8;%320000000;
sim.n_nodes = 80;
sim.n_tasks_per_node = 1; % every task needs it own copy of the fields
sim.n_cpus_per_task = 32; % these can share memory between them
sim.n_tasks = sim.n_nodes*sim.n_tasks_per_node;
sim.n_cpus = sim.n_tasks*sim.n_cpus_per_task;
sim.tot_num_particles = sim.n_tasks*sim.np_per_task;

sim.n_particle_files = sim.n_tasks;

% Memory
% Each task/proc needs the field files and the particle files
% Both fields and particle quantities are declared as 'real' in globalshf.h.
% 1 byte = 8 bits
% 1 byte = 2^30 GB ?
% 1 byte = 10^6 GB ?
sim.byte_per_GB = 10^9;2^30;

% Grid 6400x3200 gives files of 2.7 GB
% Quantities stored: (6+10*nsp)*nx*nz
% bx,by,bz,ex,ey,ez = 6
% dns,vxs,vys,vzs,pxx,pyy,pzz,pxy,pxz,pyz = 10*nsp
sim.byte_per_entry_fields = 4; % 128 bits
sim.n_entries_field_files = (6+10*sim.nsp)*sim.nx*sim.nz;
sim.size_field_files_bytes = sim.byte_per_entry_fields*sim.n_entries_field_files;
sim.size_field_files_GB = sim.size_field_files_bytes/sim.byte_per_GB;

% np = 320000000 gives files of 6.4 GB
% Quantities stored: np*5
% x,z,vx,vy,vz
sim.byte_per_entry_particles = 4; % 32 bits
sim.n_entries_particle_files = 5*sim.np_per_task;
sim.size_particle_files_bytes = sim.byte_per_entry_particles*sim.n_entries_particle_files;
sim.size_particle_files_GB = sim.size_particle_files_bytes/sim.byte_per_GB;

% Fram memory (min/max) = 64/512 (512 is for some special nodes)
sim.additional_memory_overhead_percentage = 100;
sim.memory_usage_per_node_GB = sim.n_tasks_per_node*(sim.size_particle_files_GB + sim.size_field_files_GB)*(1 + sim.additional_memory_overhead_percentage/100);
sim.memory_per_node_fram = 64;
sim.memory_per_node_vilje = 32;

% Particles per cell 
% Division of particles between species
sim.n_part(1) = sim.tot_num_particles*(0.5/9);
sim.n_part(2) = sim.tot_num_particles*(0.5/9);
sim.n_part(3) = sim.tot_num_particles*(2/9);
sim.n_part(4) = sim.tot_num_particles*(2/9);
sim.n_part(5) = sim.tot_num_particles*(2/9);
sim.n_part(6) = sim.tot_num_particles*(2/9);
sim.n_part_per_cell_average = round(sim.n_part/sim.nx/sim.nz);
% I'm not sure this meets up with the actual number of particles assigned

% More compact structure with most relevant properties
fields_to_copy = {'mime','wpewce','nx','nz','dtwpe','dtwci',...
  'dxe','dze','dxi','dzi','ximin','ximax','zimin','zimax',...
  'np_per_task','n_nodes',...
  'n_tasks_per_node','n_tasks','n_cpus_per_task','n_cpus',...
  'size_field_files_GB','size_particle_files_GB','additional_memory_overhead_percentage',...
  'memory_usage_per_node_GB',...
  'memory_per_node_fram','n_part_per_cell_average'};
for ifield = 1:numel(fields_to_copy)
  sim_compact.(fields_to_copy{ifield}) = sim.(fields_to_copy{ifield});
end

%sim
sim_compact