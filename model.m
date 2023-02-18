clc
clear

experiment = 'A';
exp_sub = 'W';
case_name = strcat(experiment, exp_sub);

Vmax = 40; U = -8; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Hmax = 3000;
dx = 5000*2; dy = 5000*2; %%%%%%%%%
nx = 501; ny = 501;
g = 9.8; f = 7.29*10^(-5);
xc = 341; yc = 251; 
%xc = 341; yc = 251 + 20;  %%%%%%%%%%%% UPPER
%xc = 341; yc = 251 - 20;  %%%%%%%%%%%% LOWER

Hm = zeros(nx); 
switch exp_sub
    case 'W'
        Sx = 50e3;  %%
        Sy = 200e3;
    case 'N'
        Sx = 200e3;
        Sy = 50e3;
    case 'S'
        R = 150e3;
    case 'L'
        R = 300e3;
end
switch experiment
    case 'A'
        Xm = 231;
        Ym = 251;
        R = 150*1e3;
        for i=1:nx   
            for j=1:nx
                Hm(j,i) = Hmax*exp(-((i-Xm)*dx)^2/Sx^2-((j-Ym)*dy)^2/Sy^2);
            end
        end 
    case 'B'
        Xm1 = 231;
        Ym1 = 201;
        Xm2 = 231;
        Ym2 = 301;
        Xm3 = 181;
        Ym3 = 251;
        Sx = 150e3;
        Sy = 150e3;
        for i=1:nx   
            for j=1:nx
                Hm(j,i) = Hmax*exp(-((i-Xm1)*dx)^2/Sx^2-((j-Ym1)*dy)^2/Sy^2) ...
                    + Hmax*exp(-((i-Xm2)*dx)^2/Sx^2-((j-Ym2)*dy)^2/Sy^2) ...
                    + Hmax*exp(-((i-Xm3)*dx)^2/Sx^2-((j-Ym3)*dy)^2/Sy^2);
            end
        end 
end

u = zeros(nx); 
v = zeros(nx);
%------tangential velocity Vt-------
gamma = 2; % Gama %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cor_xc = (xc - 1)*dx;
cor_yc = (yc - 1)*dy;
for i=1:nx 
    for j=1:ny
        cor_x = (j - 1)*dx;
        cor_y = (i - 1)*dy;
        radius = sqrt((cor_x - cor_xc)^2 + (cor_y - cor_yc)^2);
        Vt = Vmax*(radius/R)*exp(1/gamma*(1 - (radius/R)^gamma));
        vec_u = [(j-xc) (i-yc)];
        vec_vt = [-(i-yc) (j-xc)];
        vec_vt = vec_vt/norm(vec_vt);
        u(i,j) = Vt*vec_vt(1) + U;
        v(i,j) = Vt*vec_vt(2);
    end
end
u(yc,xc) = 0;
v(yc,xc) = 0;

% solve h from radial wind by eq(A4)
H = 5e3;
Ux = diff_x(u, dx);
Uy = diff_y(u, dy);
Vx = diff_x(v, dx);
Vy = diff_y(v, dy);
FG = 1/g*(f*(Vx - Uy) + 2*(Ux.*Vy - Uy.*Vx));
F = f*v/g;
G = -f*u/g;
omega = 0.95;

h = zeros(ny, nx);
h_new = h;
count = 0;

[X,Y] = meshgrid(1:10:nx,1:10:ny);
colormap(flipud(gray));
contourf(Hm);
hold on;
quiver(X, Y, u(1:10:end,1:10:end), v(1:10:end,1:10:end));
hold off;

%%

while 1 > 0
    count = count +1 ;
    disp(count)
    %implement BC at the 4 corners
    h_new(1,1) = h(1,1) + omega/4*(2*h(1,2) - 2*dx*F(1,1) + 2*h(2,1) - 2*dx*G(1,1) - 4*h(1,1)- dx^2*FG(1,1));
    h_new(1,nx) = h(1,nx) + omega/4*(2*h(1,nx-1) + 2*dx*F(1,nx) + 2*h(2,nx) - 2*dx*G(1,nx) - 4*h(1,nx) - dx^2*FG(1,nx));
    h_new(ny,1) = h(ny,1) + omega/4*(2*h(ny,2) - 2*dx*F(ny,1) + 2*h(ny-1,1) + 2*dx*G(ny,1) - 4*h(ny,1) - dx^2*FG(ny,1));
    h_new(ny,nx) = h(ny,nx) + omega/4*(2*h(ny,nx-1) + 2*dx*F(ny,nx) + 2*h(ny-1,nx) + 2*dx*G(ny,nx) - 4*h(ny,nx) - dx^2*FG(ny,nx));
    %implement BC along the 4 edges
    h_new(2:ny-1,1) = h(2:ny-1,1) + omega/4*(2*h(2:ny-1,2) - 2*dx*F(2:ny-1,1) + h(1:ny-2,1) + h(3:ny,1) -4*h(2:ny-1,1) - dx^2*FG(2:ny-1,1));
    h_new(2:ny-1,nx) = h(2:ny-1,nx) + omega/4*(2*h(2:ny-1,nx-1) + 2*dx*F(2:ny-1,nx) + h(1:ny-2,nx) + h(3:ny,nx) - 4*h(2:ny-1,nx) - dx^2*FG(2:ny-1,nx));
    h_new(1,2:nx-1) = h(1,2:nx-1) + omega/4*(h(1,1:nx-2) + h(1,3:nx) + 2*h(2,2:nx-1) - 2*dx*G(1,2:nx-1) - 4*h(1,2:nx-1) - dx^2*FG(1,2:nx-1));
    h_new(ny,2:nx-1) = h(ny,2:nx-1) + omega/4*(h(ny,1:nx-2) + h(ny,3:nx) + 2*h(ny-1,2:nx-1) + 2*dx*G(ny,2:nx-1) - 4*h(ny,2:nx-1) - dx^2*FG(ny,2:nx-1));
    %inside
    h_new(2:ny-1,2:nx-1) = h(2:ny-1,2:nx-1) + omega/4*(h(1:ny-2,2:nx-1) + h(3:ny,2:nx-1) + h(2:ny-1,1:nx-2) + h(2:ny-1,3:nx) - 4*h(2:ny-1,2:nx-1) - dx^2*FG(2:ny-1,2:nx-1));
    %compare h_new and h
    residual = abs((h(:,:) - h_new(:,:))./mean(h(:)));
    h = h_new;
    %disp(max(max(residual)));
    if residual < 1e-2 | count > 200000
        break;
    end
end
h = h + H;

case_name = 'AW_24U8';
%Hm = Hm.*0; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% No terrain 

dir = strcat('./TC_result/',case_name,'/');
dlmwrite(strcat(dir,'U/U_0.csv'), u);
dlmwrite(strcat(dir,'V/V_0.csv'), v);
dlmwrite(strcat(dir,'H/H_0.csv'), h);

disp(dir);

% Loop
dt = 24; 
 
us = zeros(nx); 
vs = zeros(nx); 
hs = zeros(nx);
un = zeros(nx); 
vn = zeros(nx); 
hn = zeros(nx);

step = 96*3600/dt;
 
for n=1:step
    us(:,:) = Update_1(u,u,dt,dx,dy);
    us = Update_2(us,v,dt,dx,dy);
    vs(:,:) = Update_1(v,u,dt,dx,dy);
    vs = Update_2(vs,v,dt,dx,dy);
    
    un(1:end,2:ny-1) = us(1:end,2:ny-1) - dt*g/(2*dx)*(h(1:end,3:ny) - h(1:end,1:ny-2)) + ...
                                                dt*f*vs(1:end,2:ny-1); 
    un(1:end,1) = un(1:end,3);
    un(1:end,ny) = un(1:end,ny-2);
    vn(2:ny-1,1:end) = vs(2:ny-1,1:end) - dt*g/(2*dy)*(h(3:ny,1:end) - h(1:ny-2,1:end)) - ...
                                                dt*f*us(2:ny-1,1:end);
    vn(1,1:end) = vn(3,1:end);
    vn(ny,1:end) = vn(ny-2,1:end);
    
    hs = Update_h1(h,h-Hm,u,dt,dx,dy);
    hs = Update_2(hs,v,dt,dx,dy);

    hn(2:ny-1,2:nx-1) = hs(2:ny-1,2:nx-1) - dt*(h(2:ny-1,2:nx-1) - Hm(2:ny-1,2:nx-1)).* ...
        ( (un(2:ny-1,3:nx) - un(2:ny-1,1:nx-2))/(2*dx) + (vn(3:ny,2:nx-1) - vn(1:ny-2,2:nx-1))/(2*dy) );
    hn(2:end-1,1) = hn(2:end-1,2);
    hn(2:end-1,end) = hn(2:end-1,end-1);
    hn(1,1:end) = hn(2,1:end);
    hn(end,1:end) = hn(end-1,1:end);
%     %----------------------------------------------------------------------
%     h = smooth(hn); 
%     u = smooth(un); 
%     v = smooth(vn);
    h = hn;
    u = un;
    v = vn;
    
    fprintf('Step: ------------------------- %d\n', n);
    
    if mod(n,3600/dt) == 0
        %fprintf('%d, v = %f, u = %f, h = %f \n', i, v_res(251,356), u_res(251,356), h_res(251,356));
        fu = [dir 'U/U_' num2str(n) '.csv'];
        fv = [dir 'V/V_' num2str(n) '.csv'];
        fh = [dir 'H/H_' num2str(n) '.csv'];

        dlmwrite(fu, u);
        dlmwrite(fv, v);
        dlmwrite(fh, h);
        
        [X,Y] = meshgrid(1:10:nx,1:10:ny);
        colormap(flipud(gray));
        contourf(Hm);
        hold on;
        quiver(X, Y, u(1:10:end,1:10:end), v(1:10:end,1:10:end));
        hold off;
        saveas(gcf,strcat('./TC_result/',case_name,'/',case_name,'_uv_', num2str(n)), 'png');
        if n ~= step
            clf
        end
    end
    
end