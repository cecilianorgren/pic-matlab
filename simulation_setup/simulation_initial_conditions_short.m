%% Constant v
nSimp = 100;
mime = 100;
Lhe = 20; % Harris sheet width in de
Lhi = Lhe/sqrt(mime); % Harris sheet width in di
l = Lhi; % harris sheet width
zh = l;
B0 = 1; % harris sheet asymptotic amplitude
BG = 0*0.25;
ah = B0*l; % rot(A) = B -> A0/l = B0 ....... rot(B) = J -> B0/l = Jy0
xp = 2.0*l; % scale length of perturbation
zp = 1.5*l; % scale length of perturbation
ap = ah*0.7; % amplitude of perturbation

TeTi = 1/5;
Ttot = 0.5;

syms x y z %ah zh ap xp zp
R = [x y z];

AH = [0; -ah*log(cosh(z/zh)); 0]; 
A = AH;
B = curl(A,R);
J = curl(B,R);
divB = divergence(B,R);
PB = 0.5*(B(1).^2+B(2).^2+B(3).^2);
% Given T, n should be calculated such that PB+PT=const
% 
PT = 0.5*B0^2-PB;
% PT = PTe + PTi = PTi/(Ti/Te) + PTi = PTi*(Te/Ti + 1)
%                = PTe + PTe/(Te/Ti) = PTe*(1 + Ti/Te)
Pi = PT/(TeTi+1);
Pe = PT/(1+1/TeTi);
n = PT/Ttot;
nH = cosh(z/zh)^-2;
% simplify(n-nH)

% Momentum conservation
%   mi*ni*vi + me*ne*ve = 0
%   ni*vi = - (me/mi)*ne*ve
% Current
%   J = ni*vi-ne*ve 
%     = ni*vi + (mi/me)*ni*vi = (1 + mi/me)*ni*vi
%     = - (me/mi)*ne*ve - ne*ve = - (me/mi + 1)*ne*ve
% Velocities
%   vi = J/[ni*(1 + mi/me)]
%   ve = -J/[ne*(me/mi + 1)]
Vi = J./(n*(mime+1));
Ve = -J./(n*(1/mime+1));
VixB = simplify(cross(Vi,B),nSimp);
VexB = cross(Ve,B);
% E+vixB - div(Pti)/ne = 0
GradPi = gradient(Pi,R);
GradPe = gradient(Pe,R);    
GradPin = GradPi/n;
GradPen = GradPe/n;
% Maybe we have to assign V from a given E (and gradPi/ne)
Ei = -VixB + GradPi./n;
Ee = -VexB - GradPe./n;

Ji = Vi*n;
Je = Ve*n;
% simplify(Vi(2))
% simplify(Ve(2))
% simplify(VixB(3))
% simplify(GradPi(3))
% simplify(n)
% simplify(GradPin(3))
% simplify(E)

f_B = matlabFunction(B(1));
f_J = matlabFunction(J(2));
f_Ji = matlabFunction(Ji(2));
f_Je = matlabFunction(Je(2));
f_PB = matlabFunction(PB);
f_PT = matlabFunction(PT);
f_Pi = matlabFunction(Pi);
f_Pe = matlabFunction(Pe);
f_Vi = matlabFunction(Vi(2));
f_Ve = matlabFunction(Ve(2));
f_VixB = matlabFunction(VixB(3));
f_VexB = matlabFunction(VexB(3));
f_GradPi = matlabFunction(GradPi(3));
f_GradPe = matlabFunction(GradPe(3));
f_GradPin = matlabFunction(GradPin(3));
f_GradPen = matlabFunction(GradPen(3));
f_Ei = matlabFunction(Ei(3));
f_Ee = matlabFunction(Ee(3));
f_N = matlabFunction(n);

% Plot
zvec = linspace(-5,5,100);
nrows = 8; ncols = 1; npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;
if 1 % Bx, Jy, n
  hca = h(isub); isub = isub + 1;  
  plot(hca,zvec,f_B(zvec),zvec,f_J(zvec),zvec,f_N(zvec)) 
  legend(hca,{'B_x','J_y','n'})
end
if 1 % P
  hca = h(isub); isub = isub + 1;  
  plot(hca,zvec,f_PB(zvec),zvec,f_PT(zvec),zvec,f_Pi(zvec),zvec,f_Pe(zvec),zvec,f_PB(zvec)+f_PT(zvec))  
  legend(hca,{'P_B','P_T=P_i+P_e','P_i','P_e','P_B+P_T'})
end
if 1 % V
  hca = h(isub); isub = isub + 1;  
  plot(hca,zvec,f_Vi(zvec),zvec,f_Ve(zvec))  
  legend(hca,{'v_i','v_e'})
  hca.YLim = hca.YLim+diff(hca.YLim)*0.05*[-1 1];
end
if 1 % Je,Ji (flux)
  hca = h(isub); isub = isub + 1;  
  plot(hca,zvec,f_Ji(zvec),zvec,f_Je(zvec))  
  legend(hca,{'j_i','j_e'})
  hca.YLim = hca.YLim+diff(hca.YLim)*0.05*[-1 1];
end
if 1 % VxB
  hca = h(isub); isub = isub + 1;  
  plot(hca,zvec,f_VixB(zvec),zvec,f_VexB(zvec))  
  legend(hca,{'v_ixB','v_exB'})  
end
if 1 % GradP
  hca = h(isub); isub = isub + 1;  
  plot(hca,zvec,f_GradPi(zvec),zvec,f_GradPe(zvec))  
  legend(hca,{'grad(P_i)','grad(P_e)'})  
end
if 1 % GradP/n
  hca = h(isub); isub = isub + 1;  
  plot(hca,zvec,f_GradPin(zvec),zvec,f_GradPen(zvec))  
  legend(hca,{'grad(P_i)/n','grad(P_e)/n'})  
end
if 1 % E
  hca = h(isub); isub = isub + 1;  
  plot(hca,zvec,f_Ei(zvec),zvec,f_Ee(zvec),'--')  
  legend(hca,{'-v_ixB+grad(P_i)/n','-v_exB-grad(P_e)/n'})  
  hca.YLabel.String = 'E_z';
end
compact_panels(0.01)
for ip = 1:npanels
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
end

%% Get v's from gradP/n and E
nSimp = 100;
mime = 100;
Lhe = 20; % Harris sheet width in de
Lhi = Lhe/sqrt(mime); % Harris sheet width in di
l = Lhi; % harris sheet width
zh = l;
B0 = 1; % harris sheet asymptotic amplitude
BG = 0*0.25;
ah = B0*l; % rot(A) = B -> A0/l = B0 ....... rot(B) = J -> B0/l = Jy0
xp = 2.0*l; % scale length of perturbation
zp = 1.5*l; % scale length of perturbation
ap = ah*0.7; % amplitude of perturbation

TeTi = 1/5;
Ttot = 0.5;

syms x y z %ah zh ap xp zp
R = [x y z];

AH = [0; -ah*log(cosh(z/zh)); 0]; 
A = AH;
B = curl(A,R);
J = curl(B,R);
divB = divergence(B,R);
PB = 0.5*(B(1).^2+B(2).^2+B(3).^2);
% Given T, n should be calculated such that PB+PT=const
% 
PT = 0.5*B0^2-PB;
% PT = PTe + PTi = PTi/(Ti/Te) + PTi = PTi*(Te/Ti + 1)
%                = PTe + PTe/(Te/Ti) = PTe*(1 + Ti/Te)
Pi = PT/(TeTi+1);
Pe = PT/(1+1/TeTi);
n = PT/Ttot;
nH = cosh(z/zh)^-2;
phi = -0.5*n;
phi = -0.5*cosh(z/(zh*1))^-2;
E = -gradient(phi,z);
% simplify(n-nH)

% Momentum conservation
%   mi*ni*vi + me*ne*ve = 0
%   ni*vi = - (me/mi)*ne*ve
% Current
%   J = ni*vi-ne*ve 
%     = ni*vi + (mi/me)*ni*vi = (1 + mi/me)*ni*vi
%     = - (me/mi)*ne*ve - ne*ve = - (me/mi + 1)*ne*ve
% Velocities
%   vi = J/[ni*(1 + mi/me)]
%   ve = -J/[ne*(me/mi + 1)]
% Vi = J./(n*(mime+1));
% Ve = -J./(n*(1/mime+1));
% VixB = simplify(cross(Vi,B),nSimp);
% VexB = cross(Ve,B);
% E+vixB - div(Pti)/ne = 0
GradPi = gradient(Pi,R);
GradPe = gradient(Pe,R);    
GradPin = GradPi/n;
GradPen = GradPe/n;
% Maybe we have to assign V from a given E (and gradPi/ne)
VixB = -E + GradPi./n; % E+vixB-gradPi/n
VexB = -E - GradPe./n; % -E-vexB-gradPe/n
Vi = -VixB(3)/B(1); % vB_z = (vxBy-vyBx) ~ +vyBx -> vy = -vB_z/Bx
Ve = -VexB(3)/B(1); %

Ji = Vi*n;
Je = Ve*n;
% simplify(Vi(2))
% simplify(Ve(2))
% simplify(VixB(3))
% simplify(GradPi(3))
% simplify(n)               
% simplify(GradPin(3))
% simplify(E)

f_Phi = matlabFunction(phi);
f_B = matlabFunction(B(1));
f_J = matlabFunction(J(2));
f_Ji = matlabFunction(Ji);
f_Je = matlabFunction(Je);
f_PB = matlabFunction(PB);
f_PT = matlabFunction(PT);
f_Pi = matlabFunction(Pi);
f_Pe = matlabFunction(Pe);
f_Vi = matlabFunction(Vi);
f_Ve = matlabFunction(Ve);
f_VixB = matlabFunction(VixB(3));
f_VexB = matlabFunction(VexB(3));
f_GradPi = matlabFunction(GradPi(3));
f_GradPe = matlabFunction(GradPe(3));
f_GradPin = matlabFunction(GradPin(3));
f_GradPen = matlabFunction(GradPen(3));
f_E = matlabFunction(E);
f_N = matlabFunction(n);

int_Ji = trapz(zvec,f_Ji(zvec));
int_Je = trapz(zvec,f_Je(zvec));
int_Flux = int_Ji + int_Je;

% Plot
zvec = linspace(-10,10,100);
nrows = 8; ncols = 1; npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;
if 1 % Bx, Jy, n
  hca = h(isub); isub = isub + 1;  
  plot(hca,zvec,f_B(zvec),zvec,f_J(zvec),zvec,f_N(zvec)) 
  legend(hca,{'B_x','J_y','n'})
end
if 1 % P
  hca = h(isub); isub = isub + 1;  
  plot(hca,zvec,f_PB(zvec),zvec,f_PT(zvec),zvec,f_Pi(zvec),zvec,f_Pe(zvec),zvec,f_PB(zvec)+f_PT(zvec))  
  legend(hca,{'P_B','P_T=P_i+P_e','P_i','P_e','P_B+P_T'})
end
if 0 % phi
  hca = h(isub); isub = isub + 1;  
  plot(hca,zvec,f_Phi(zvec)) 
  legend(hca,{['\phi = ' char(phi)]})
end
if 1 % phi, E
  hca = h(isub); isub = isub + 1;  
  plot(hca,zvec,f_Phi(zvec),zvec,f_E(zvec)) 
  legend(hca,{['\phi = ' char(phi)],['E = -\partial_z\phi = ' char(E)]},'Box','off','Location','best')
end
if 1 % GradP
  hca = h(isub); isub = isub + 1;  
  plot(hca,zvec,f_GradPi(zvec),zvec,f_GradPe(zvec))  
  legend(hca,{'grad(P_i)','grad(P_e)'})  
end
if 1 % GradP/n
  hca = h(isub); isub = isub + 1;  
  plot(hca,zvec,f_GradPin(zvec),zvec,f_GradPen(zvec))  
  legend(hca,{'grad(P_i)/n','grad(P_e)/n'})  
end
if 1 % VxB
  hca = h(isub); isub = isub + 1;  
  plot(hca,zvec,f_VixB(zvec),zvec,f_VexB(zvec))  
  legend(hca,{'v_ixB','v_exB'})  
end
if 1 % V
  hca = h(isub); isub = isub + 1;  
  plot(hca,zvec,f_Vi(zvec),zvec,f_Ve(zvec))  
  legend(hca,{'v_i','v_e'})
  hca.YLim = hca.YLim+diff(hca.YLim)*0.05*[-1 1];
end
if 1 % Je,Ji (flux)
  hca = h(isub); isub = isub + 1;  
  plot(hca,zvec,f_Ji(zvec),zvec,f_Je(zvec))  
  legend(hca,{sprintf('j_i: int j_i dz = %.2f',int_Ji),sprintf('j_e: int j_e dz = %.2f',int_Je)},'Box','off','Location','best')
  hca.YLim = hca.YLim+diff(hca.YLim)*0.05*[-1 1];
end
if 0 % E
  hca = h(isub); isub = isub + 1;  
  plot(hca,zvec,f_E(zvec))  
  legend(hca,{'-grad(phi)'})  
  hca.YLabel.String = 'E_z';
end
compact_panels(0.01)
for ip = 1:npanels
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
end