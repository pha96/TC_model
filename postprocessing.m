clc
clear all;

experiment = 'A';
exp_sub = 'W';
case_name = strcat(experiment, exp_sub);
Hmax = 3000;
dx = 5000*2; dy = 5000*2;
nx = 501; ny = 501;
f = 7.29*10^(-5);

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
nstep = 96*3600/12; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
step = nstep/600 + 1;

%%

%case_name = 'BL_8';
xc = zeros(step,1);
yc = zeros(step,1);

dir = strcat('./TC_result/',case_name,'/');

vor = zeros(step,nx,ny);
st = zeros(step,nx,ny);
vmax = zeros(step,1);
u = vor;
v = u;
pv = v;

dt = 12;
for n = 1:step
    disp(n);
    fu = [dir 'U/U_' num2str((n-1)*3600/dt) '.csv'];
    fv = [dir 'V/V_' num2str((n-1)*3600/dt) '.csv'];
    fh = [dir 'H/H_' num2str((n-1)*3600/dt) '.csv'];
    
    u(n,:,:) = csvread(fu);
    v(n,:,:) = csvread(fv);
    h = csvread(fh);
    
    uu = csvread(fu);
    vv = csvread(fv);
    
    [yy, xx] = find(h == min(h(:)));
    yc(n) = yy(1);
    xc(n) = xx(1);
    
    uv = sqrt(u(n,:,:).^2 + v(n,:,:).^2);
    vmax(n) = max(uv(:));
    
%     if (mod(n,6) == 1)
%         vr = diff_x(vv,dx) - diff_y(uu,dy);
%         vor(n,:,:) = vr;
%         st(n,:,:) = stream_function(vv, -uu, vr, dx, dy);
%         pv(n,:,:) = (vr + f)./h;
%         csvwrite([dir 'ST_' num2str(n) '.csv'], st(n,:,:));
%         csvwrite([dir 'PV_' num2str(n) '.csv'], pv(n,:,:));
%     end
end

%% Storm track
colormap(flipud(gray));
contourf(Hm);
hold on;
[X,Y] = meshgrid(1:10:nx,1:10:ny);
quiver(X, Y, uu(1:10:end,1:10:end), vv(1:10:end,1:10:end),'b');
plot(xc(2:end), yc(2:end),'r','LineWidth',1.5);
file_name = strcat('Storm track',{' '},case_name);
title(case_name);
xlim([0 501]);
ylim([0 501]);
saveas(gcf,char(strrep(file_name,' ','')), 'png');

%% Vmax
plot(vmax, 'LineWidth', 1.5);
grid on;
xlabel('Hours');
xlim([0 96]);
%ylim([0 50]);
ylabel('Vmax (m/s)');
title(case_name);
file_name = strcat('Vmax','_',case_name);
saveas(gcf,char(file_name), 'png');

%% Asymmetric wind
for n = 7:6:step
    i = fix((n-7)/6) + 1;
    subplot(4,4,i);
    center = yc(n);
    left = xc(n) - 99;
    right = xc(n) + 99;
    uv = sqrt(u(n,:,:).^2 + v(n,:,:).^2);
    plot(squeeze(uv(1,center,left:right)),'LineWidth',1.5);
    xticks([0 50 100 150 200]);
    xticklabels({'-1000','-500','0','50','1000'});
    if (mod(i,4) == 1)
        ylabel('Wind speed (m/s)','FontWeight','bold');
    end
    if (i >= 13 & i <= 16)
        xlabel('Distance (km)', 'FontWeight','bold');
    end
    xlim([0 200]);
    ylim([0 80]);
    title([num2str(n-1) ' hour']);
    grid on;
end
file_name = strcat('CrossWind_',case_name);
set(gcf,'units','normalized','outerposition',[0 0 1 1])
saveas(gcf,char(file_name), 'png');

%% Vector
for n = 7:6:step
    i = fix((n-7)/6) + 1;
    subplot(4,4,i);
    [X,Y] = meshgrid(1:20:nx,1:20:ny);
    hold on;
    quiver(X, Y, squeeze(u(n,1:20:end,1:20:end)), squeeze(v(n,1:20:end,1:20:end)),'b');
    contour(Hm, 'k');
    title([num2str(n-1) ' hour']);
    xlim([0 501]);
    ylim([0 501]);
    hold off;
end
file_name = strcat('Wind_',case_name);
set(gcf,'units','normalized','outerposition',[0 0 1 1])
saveas(gcf,char(file_name), 'png');

%% Stream Function
for n = 7:6:step
    i = fix((n-7)/6) + 1;
    subplot(4,4,i);
    contourf(squeeze(st(n,:,:)), 'edgecolor','none');
    hold on;
    contour(Hm, 'k');
    title([num2str(n-1) ' hour']);
    xlim([0 501]);
    ylim([0 501]);
    colorbar
    hold off;
end
file_name = strcat('ST_',case_name);
set(gcf,'units','normalized','outerposition',[0 0 1 1])
saveas(gcf,char(file_name), 'png');

%% Potential Velocity
for n = 7:6:step
    i = fix((n-7)/6) + 1;
    subplot(4,4,i);
    contourf(squeeze(pv(n,:,:))*1e8, 'edgecolor','none');
    hold on;
    contour(Hm, 'k');
    title([num2str(n-1) ' hour']);
    xlim([0 501]);
    ylim([0 501]);
    colorbar;
    grid on;
    hold off;
end
file_name = strcat('PV_',case_name);
set(gcf,'units','normalized','outerposition',[0 0 1 1])
saveas(gcf,char(file_name), 'png');

%% Sensitivity track
exp = ["AW_g05" "AW_g1" "AW"];
track_x = zeros(3,step);
track_y = zeros(3,step);

hold on;
[X,Y] = meshgrid(1:20:nx,1:20:ny);
colormap(flipud(gray));
contourf(Hm);

plt = zeros(1,3);
for i = 1:length(plt)
    dir = strcat('./TC_result/',exp(i),'/');
    disp(i);
    for n = 1:step
        fh = char(strcat(dir,'H/H_',num2str((n-1)*600),'.csv'));
        h = csvread(fh);
        
        [yy, xx] = find(h == min(h(:)));
        track_y(i,n) = yy(1);
        track_x(i,n) = xx(1);
    end
    plt(i) = plot(squeeze(track_x(i,:)), squeeze(track_y(i,:)),'LineWidth',1.5);
    %hold on;
end

file_name = 'Track_gamma';
title('Sensitivity Track');
xlim([0 501]);
ylim([0 501]);
legend(plt(:), 'G0.5', 'G1', 'G2');
saveas(gcf,file_name, 'png');

%% Sensitivity Vmax
exp = ["AW_g05" "AW_g1" "AW"];
v_exp = zeros(3,step);

hold on;
for i = 1:length(exp)
    dir = strcat('./TC_result/',exp(i),'/');
    disp(i);
    for n = 1:step
        fu = char(strcat(dir,'U/U_',num2str((n-1)*600),'.csv'));
        fv = char(strcat(dir,'V/V_',num2str((n-1)*600),'.csv'));

        uu = csvread(fu);
        vv = csvread(fv);     
        uv = sqrt(uu.^2 + vv.^2);
        v_exp(i,n) = max(uv(:));      
    end
    plot(squeeze(v_exp(i,:)), 'LineWidth', 1.5);
end

file_name = 'Vmax_exp';
grid on;
xlabel('Hours', 'FontWeight','bold');
xlim([0 96]);
ylabel('Vmax (m/s)', 'FontWeight','bold');
legend('G0.5', 'G1', 'G2', 'Location','northwest');
title('Sensitivity Vmax');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
saveas(gcf,file_name, 'png');

%% Wind speed zooming
for n = 55:6:step
    i = fix((n-55)/6) + 1;
    subplot(2,4,i);
    uv = sqrt(u(n,:,:).^2 + v(n,:,:).^2);
    contourf(squeeze(uv), 'edgecolor','none');
    colormap(jet);
    caxis([0 70]);
    [X,Y] = meshgrid(150:10:350,150:10:350);
    hold on;
    quiver(X, Y, squeeze(u(n,150:10:350,150:10:350)), squeeze(v(n,150:10:350,150:10:350)),'k');
    contour(Hm, 'k');
    title([num2str(n-1) ' hour']);
    xlim([150 350]);
    ylim([150 350]);
    hold off;
end
file_name = strcat('W_',case_name);
set(gcf,'units','normalized','outerposition',[0 0 1 1])
saveas(gcf,char(file_name), 'png');