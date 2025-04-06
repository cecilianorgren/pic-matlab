function  x_res = eom_pic_external(t,xyx_vxvyvz,fEx,fEy,fEz,fBx,fBy,fBz,m,q)
    x = xyx_vxvyvz(1);
    y = xyx_vxvyvz(2);
    z = xyx_vxvyvz(3);
    vx = xyx_vxvyvz(4);
    vy = xyx_vxvyvz(5);
    vz = xyx_vxvyvz(6);

    if isnan(x)
      1;
    end
    %disp(sprintf('t = %g, x = %g, y = %g, z = %g',t,x,y,z))

    %method = 'spline';        

    %[Ex,Ey,Ez,Bx,By,Bz] = pic.interp_EB(x,z,t,method);
    Ex = fEx(x,z);
    Ey = fEy(x,z);
    Ez = fEz(x,z);
    Bx = fBx(x,z);
    By = fBy(x,z);
    Bz = fBz(x,z);
    
    %disp(sprintf('%.3f %.3f %.3f, %.3f %.3f %.3f, %.3f %.3f %.3f, %.3f, %.3f %.3f',Ex,Ey,Ez,Bx,By,Bz,x_vect(1),x_vect(2),x_vect(3),x_vect(4),x_vect(5),x_vect(6)))        
    %it = size(x_sol_all,1);
    %x_sol_all(it+1,1) = Ex;
    %x_sol_all(it+1,2) = Ey;
    %x_sol_all(it+1,3) = Ez;
    %x_sol_all(it+1,4) = Bx;
    %x_sol_all(it+1,5) = By;
    %x_sol_all(it+1,6) = Bz;
    %                                                                                                                                                 x_sol_all(it+1,7) = t;

    % Equations to be solved
    x_res = zeros(6,1);
    x_res(1) = vx; % dx/dt = vx;
    x_res(2) = vy; % dy/dt = vy;
    x_res(3) = vz; % dz/dt = vz;
    x_res(4) = (q/m)*(Ex + vy*Bz - vz*By);
    x_res(5) = (q/m)*(Ey + vz*Bx - vx*Bz);
    x_res(6) = (q/m)*(Ez + vx*By - vy*Bx);                                              

  end      