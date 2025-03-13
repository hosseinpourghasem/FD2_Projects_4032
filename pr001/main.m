clc; clear; close all; 
dbstop error

% تعریف متغیرهای سراسری
global m g Ixx Iyy Izz Ixz
global FAX FAY FAZ FTX FTY FTZ La Ma Na Lt Mt Nt

% بارگذاری داده‌ها از فایل
load data.txt

%% **تعریف تابع معادلات حرکت**
fun = @(t, states) EOM_6DoF_Body(t, states);

%% **پارامترهای جرم و اینرسی هواپیما**
Weight = 564000; % lbs
g = 32.2; % ft/s²
m = Weight / g; % محاسبه جرم
Ixx = 13700000; % slug.ft²
Iyy = 30500000; % slug.ft²
Izz = 49700000; % slug.ft²
Ixz = 830000; % slug.ft²

%% **شرایط پروازی اولیه**
Altitude = 1000; % ft
U0 = 262.6411; % ft/s (سرعت اولیه)
V0 = 0; % Symmetric flight (زاویه جانبی صفر)
W0 = 22.2674; 
P0 = 0; Q0 = 0; R0 = 0; % نرخ‌های زاویه‌ای اولیه
phi0 = 0; psi0 = 0; % زوایای اولیه
theta0 = deg2rad(1.81); % تبدیل درجه به رادیان
X0 = 0; Y0 = 0; Z0 = -Altitude; % موقعیت اولیه

% مقداردهی اولیه متغیرهای وضعیت
IC = [U0 V0 W0 P0 Q0 R0 phi0 theta0 psi0 X0 Y0 Z0]; 
states(1, :) = IC;

%% **شبیه‌سازی با روش اویلر**
t = data(:, 1); 
dt = t(2) - t(1);

for k = 2:numel(t)
    disp(['t = ' num2str(t(k))])

    % استخراج داده‌های نیروها و گشتاورها از data.txt
    FAX = data(k, 2); FAY = data(k, 3); FAZ = data(k, 4);
    FTX = data(k, 5); FTY = data(k, 6); FTZ = data(k, 7);
    La = data(k, 8);  Ma = data(k, 9);  Na = data(k, 10);
    Lt = data(k, 11); Mt = data(k, 12); Nt = data(k, 13);

    % حل معادلات با روش اویلر 
    states(k, :) = states(k-1, :) + dt * fun(t(k-1), states(k-1, :) ); 
end

%% **استخراج متغیرها**
U = states(:, 1); V = states(:, 2); W = states(:, 3);
P = states(:, 4); Q = states(:, 5); R = states(:, 6);
phi = states(:, 7); theta = states(:, 8); psi = states(:, 9);
X = states(:, 10); Y = states(:, 11); Z = states(:, 12);

%% **رسم مسیر پرواز سه‌بعدی**
fig1=figure(1);
set(fig1,'name','UVW')
subplot(411)
plot(t,U,'LineWidth',1)
grid on
xlabel('Time (sec)')
ylabel('U(ft/s)')
subplot(412)
plot(t,V,'LineWidth',1)
grid on
xlabel('Time (sec)')
ylabel('V(ft/s)')
subplot(413)
plot(t,W,'LineWidth',1)
grid on
xlabel('Time (sec)')
ylabel('W(ft/s)')
subplot(414)
plot(t,sqrt(U.^2+V.^2+W.^2),'LineWidth',1)
grid on
xlabel('Time (sec)')
ylabel('Velocity(ft/s)')

fig2=figure(2);
set(fig2,'name','PQR')
subplot(311)
plot(t,P*180/pi,'LineWidth',1)
grid on
xlabel('Time(sec)')
ylabel('P(deg/s')
subplot(312)
plot(t,Q*180/pi,'LineWidth',1)
grid on
xlabel('Time(sec)')
ylabel('Q(deg/s')
subplot(313)
plot(t,R*180/pi,'LineWidth',1)
grid on
xlabel('Time(sec)')
ylabel('R(deg/s')

fig3=figure(3);
set(fig3,'name','Euler Angles')
subplot(311)
plot(t,phi*180/pi,'LineWidth',1)
grid on
xlabel('Time (sec)')
ylabel('\phi (deg)')
subplot(312)
plot(t,theta*180/pi,'LineWidth',1)
grid on
xlabel('Time (sec)')
ylabel('\theta (deg)')
subplot(313)
plot(t,psi*180/pi,'LineWidth',1)
grid on
xlabel('Time (sec)')
ylabel('\psi (deg)')


fig4=figure(4);
set(fig4,'name','position')
subplot(331)
plot(t,X,'LineWidth',1)
grid on
xlabel('Time (sec)')
ylabel('X (ft)')
subplot(334)
plot(t,Y,'LineWidth',1)
grid on
xlabel('Time (sec)')
ylabel('Y (ft)')
subplot(337)
plot(t,Z,'LineWidth',1)
grid on
xlabel('Time (sec)')
ylabel('Z (ft)')
subplot(3,3,[2 3 5 6 8 9])
plot3(X,Y,Z,'LineWidth',1)
grid on
xlabel('X (ft)')
ylabel('Y (ft)')
zlabel('Z (ft)')


fig5=figure(5);
set(fig5,'name','Euler Angles')
subplot(211)
plot(t,atan((W./U))*180/pi,'LineWidth',1)
grid on
xlabel('Time (sec)')
ylabel('\alpha (deg)')
subplot(212)
plot(t,asin((V./sqrt(U.^2+V.^2+W.^2)))*180/pi,'LineWidth',1)
grid on 
xlabel('Time (sec)')
ylabel('\beta (deg)')

figHandles = findall(0, 'Type', 'figure'); % گرفتن همه فیگورها
filename = 'All_Plots.pdf'; % نام فایل خروجی

% چاپ همه فیگورها در یک فایل PDF
for k = 1:length(figHandles)
    figure(figHandles(k));
    if k == 1
        exportgraphics(gcf, filename, 'ContentType', 'vector', 'Append', false);
    else
        exportgraphics(gcf, filename, 'ContentType', 'vector', 'Append', true);
    end
end




function out = EOM_6DoF_Body(t, states)
    global m g Ixx Iyy Izz Ixz
    global FAX FAY FAZ FTX FTY FTZ La Ma Na Lt Mt Nt

    % استخراج وضعیت‌های فعلی
    U = states(1); V = states(2); W = states(3);
    P = states(4); Q = states(5); R = states(6);
    phi = states(7); theta = states(8); psi = states(9);
    X=states(10); Y=states(11);Z=states(12);
    % نیروهای کل وارد بر بدنه
    F_x = FAX + FTX;
    F_y = FAY + FTY;
    F_z = FAZ + FTZ;
    
    % معادلات حرکت خطی
    Udot = (1/m) * (-m*g*sin(theta) + FAX+FTX) - Q*W + R*V;
    Vdot = (1/m) * (m*g*cos(theta)*sin(phi) + FAY+FTY) - R*U + P*W;
    Wdot = (1/m) * (m*g*cos(theta)*cos(phi) + FAZ+FTZ) - P*V + U*Q;

    % معادلات حرکت زاویه‌ای
    L = La + Lt;
    M = Ma + Mt;
    N = Na + Nt;
    
    Pdot = (1 / (Ixx*Izz - Ixz^2)) * ((Ixx - Iyy + Izz) * Ixz * P * Q + (Iyy * Izz - Izz^2 - Ixz^2) * Q * R + Izz * L + Ixz * N);
    Qdot = (1 / Iyy) * ((Izz - Ixx) * P * R + (R^2 - P^2) * Ixz + M);
    Rdot = (1 / (Ixx * Izz - Ixz^2)) * ((Ixx^2 - Ixx * Iyy + Ixz^2) * P * Q + (-Ixx + Iyy - Izz) * Ixz * Q * R + Ixz * L + Ixx * N);

    % نرخ تغییر زوایا
    phidot = P + Q*sin(phi)*tan(theta) + R*cos(phi)*tan(theta);
    thetadot = Q*cos(phi) - R*sin(phi);
    psidot = (Q*sin(phi) + R*cos(phi)) / cos(theta);

    % ماتریس تبدیل از دستگاه بدنی به اینرسی
    Cb2i = [cos(psi)*cos(theta),cos(psi)*sin(phi)*sin(theta)-cos(phi)*sin(psi), sin(phi)*sin(psi)+cos(phi)*cos(psi)*sin(theta);cos(theta)*sin(psi),cos(phi)*cos(psi)+sin(phi)*sin(psi)*sin(theta),cos(phi)*sin(psi)*sin(theta)-cos(psi)*sin(phi);-sin(theta),cos(theta)*sin(phi),cos(theta)*cos(phi)];
    
    XYZdot = Cb2i * [U; V; W]; 
    Xdot=XYZdot(1);
    Ydot=XYZdot(2);
    Zdot=XYZdot(3);

    out = [Udot Vdot Wdot Pdot Qdot Rdot phidot thetadot psidot Xdot Ydot Zdot];

end