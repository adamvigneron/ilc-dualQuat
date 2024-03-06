%% REVIEWED LITERATURE

% "Iterative Learning Control of Spacecraft Proximity Operations Based on Confidence Level"
% Steve Ulrich and Kirk Hovell
% AIAA GNC 2017
% https://arc.aiaa.org/doi/10.2514/6.2017-1046 

% This script re-creates the numerical example of Section III.C


%% SCRIPT CONFIGURATION

% switch confidence level factor
% 1 - constant beta
% 2 - varying beta
clfSwitch = 2;

% choose reference trajectory
% 1 - v0 != 0 (baseline)
% 2 - v0 == 0 (proposal) 
refTraj = 2;


%% PROLOGUE

% dock all figures by default
set(0,'DefaultFigureWindowStyle','docked');

% constants
mu  = 398600;  % km^3/s^2
O33 = zeros(3,3);
I33 = eye(3);

% parameters
m    = 20;  % kg
a    = 7200;  % km
f    = 240;  % Hz
dt   = 1/f;  % s
Fmax = 2.5;  % N


%% DEFINE A REFERENCE TRAJECTORY

switch refTraj

    case 1  % baseline trajectory with v0 != 0

        tf = 300;  % s
        t  = 0:dt:tf;
        r  = 50;
        
        yd = [0*t; -r*cos(2*pi*t/tf); r*sin(2*pi*t/tf)];
        x0 = [  0;                -r;                 0; ... 
                0;                 0;                 0];
    
    case 2  % proposed trajectory with v0 == 0

        tf = 300;  % s
        r  = 50;   % m
        Fa = 1.0;  % N, allocated force

        v1 = 2*pi*r/tf;
        t1 = m*v1/Fa;

        t  = 0:dt:(300+floor(t1));

        yd = zeros(3,length(t));
        x0 = [  0;                -r;                 0; ...
                0;                 0;                 0];

        for k = 1:length(t)
            tk = t(k);

            if tk < t1
                s = 0.5 * Fa/m * tk^2;
            else
                s = 0.5 * Fa/m * t1^2 + v1 * (tk-t1);
            end

            alpha = s/r;

            yd(2,k) = -50*cos(alpha);
            yd(3,k) =  50*sin(alpha);
        end

    otherwise
        error('Invalid refTraj value');

end

% plot the reference trajectory
fig1 = figure('Name','trajectory');
plot( yd(2,:), yd(3,:), 'k--', 'LineWidth', 2 );
title('trajectory');
xlabel('y, m');
ylabel('z, m');
legend('reference');
grid on;
hold on;


%% DEFINE STATE-SPACE MODELS

% plant -- double-integrator
A = [O33 I33; ...
     O33 O33];
B = [O33; ...
     I33/m];
C = [I33 O33];
D = O33;
sys2i = c2d(ss(A,B,C,D),dt);

% plant -- clohessy-wiltshire
w = sqrt(mu/a^3);
Acw = [O33 I33; ...
       3*w^2 0    0    0 2*w 0; ...
       0     0    0 -2*w   0 0; ...
       0     0 -w^2    0   0 0];
sysCW = c2d(ss(Acw,B,C,D),dt);


% controller -- PD
kp = [0.172     0     0; ...
      0     0.172     0; ...
      0         0 0.100];

kd = [4 0 0; ...
      0 4 0; ...
      0 0 3];


%% INITIALIZE MAIN LOOP AND OUTPUT FIGURES

% loop initialization
u      = zeros(3,length(t));
x      = zeros(6,length(t));
x(:,1) = x0;
N      = 10;

% y-axis control signal figure
fig2 = figure('Name','y-control');
plot(t,u(2,:));
title('control')
xlabel('time, s');
ylabel('f_y, N');
grid on;
hold on;

% z-axis control signal figure
fig3 = figure('Name','z-control');
plot(t,u(3,:));
title('control')
xlabel('time, s');
ylabel('f_z, N');
grid on;
hold on;


%% DETAILED EXAMINATION OF Z-CONTROL

u0 = zeros(1,length(t));
uP = zeros(1,length(t));
uD = zeros(1,length(t));

fig4 = figure('Name','z-detail');

subplot(1,3,1);
plot(t,u(3,:));
title('learning')
xlabel('time, s')
ylabel('f_z, N');
grid on;
hold on;

subplot(1,3,2);
plot(t,u(3,:));
title('proportional')
xlabel('time, s')
grid on;
hold on;

subplot(1,3,3);
plot(t,u(3,:));
title('derivative')
xlabel('time, s');
grid on;
hold on;


%% MAIN LOOP

for j = 1:N  % iteration loop

    % fixed or variable confidence level factor
    switch (clfSwitch)
        case 1
            beta = 1;
        case 2
            beta = (j-1)/N;  % we use beta_{j}, not beta_{j+1}
        otherwise
            error('Invalid clfSwitch value');
    end

    for k = 1:(length(t)-1)  % timestep loop

        % control error; state was calculated in the _previous_ timestep
        e = yd(:,k) - sysCW.c * x(:,k);

        % extract data for additional analysis of z-control
        u0(k) = beta * u(3,k);
        uP(k) = kp(3,3) * e(3);
        
        % control command has three parts: learning, proportional...
        u(:,k) = beta * u(:,k) + kp * e;
 
        % ...and derivative
        if k > 1
            e0     = yd(:,k-1) - sysCW.c * x(:,k-1);

            % extract data for additional analysis of z-control
            uD(k)  = kd(3,3) * (e(3)-e0(3)) / dt;

            u(:,k) = u(:,k) + kd * (e-e0) / dt;
        end

        % apply saturation condition
        for i = 1:3  % axis loop
            u(i,k) = min( Fmax, max( -Fmax, u(i,k) ) );
        end

        % state at the _next_ timestep
        x(:,k+1) = sysCW.a * x(:,k) + sysCW.b * u(:,k);

    end

    % print average 3D position error
    posErr3D = mean( vecnorm( x(1:3,:)-yd ) );
    disp( ['Iteration ' num2str(j,'%02i') ': posErr3D = ' num2str(posErr3D) ] );

    % plot trajectory on existing figure
    set(groot,'CurrentFigure',fig1);
    plot( x(2,:), x(3,:) );

    % plot y-control on existing figure
    set(groot,'CurrentFigure',fig2);
    plot( t, u(2,:));

    % plot z-control on existing figure
    set(groot,'CurrentFigure',fig3);
    plot( t, u(3,:));

    % plot z-details on existing figure
    set(groot,'CurrentFigure',fig4);
    subplot(1,3,1); plot(t,u0);
    subplot(1,3,2); plot(t,uP);
    subplot(1,3,3); plot(t,uD);

end


%% POSTSCRIPT

% reset figure style
set(0,'DefaultFigureWindowStyle','normal')

