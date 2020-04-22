%illustrate_timing_limits

% Two spacecraft located at 
x1 = 10; % km
x2 = 20; % km


% instrument sampling frequency
f = 100; % Hz, s^-1

% Structure
% speed of structure 
v_str = 100; % km/s
% "size" of structure
l_str = 100; % km
% location of structure at t = 0;
x0_str = 3*l_str;
% structure of structure
data = @(x_sc,t,v_str,x0_str) tanh((t*v_str-x0_str-x_sc)/l_str);

data1 = data(x1,t,v_str,x0_str);
data2 = data(x2,t,v_str,x0_str);

% cross correlation to get the time delay
[acor,lag] = xcorr(data1,data2,'coeff');

nstop = 2*x0_str/v_str; % s
t = linspace(0,nstop,nstop*f);

hca = subplot(2,1,1);
plot(hca,t,data1,'.-',t,data2,'.-')

hca = subplot(2,1,2);
plot(hca,lag,acor)
