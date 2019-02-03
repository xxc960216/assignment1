function [] = Q1(numOfAtom, numOfStep)
% ELEC4700 - Assignment 1
% Xiaochen Xin 100989338
C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass(kg)
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665; %metres (32.1740 ft) per s²
C.am = 1.66053892e-27;

%1-1
mn = 0.26*C.m_0; %Effective Mass
l = 200e-9; %Length of area (m)
w = 100e-9; %Width of area (m)
T = 300; %Kelvin

vth = sqrt(C.kb*T/mn)%thermal velocity(velocity at which the particles are travelling at)

%1-2
tmn = 0.2e-12; %mean time between collision (s)
mfp = tmn* vth
 
%1-3
xr = 200e-9.*rand(numOfAtom,1); %x of 100 random locations
yr = 100e-9.*rand(numOfAtom,1); %y of 100 random locations
%Define two arrays store the previous locations
xrp = xr;
yrp = yr;
ang = 2*pi.*rand(numOfAtom,1); %angle in rad of 100 random locations
vx = cos(ang)*vth; %initial horizontal velocity
vy = sin(ang)*vth; %initial vertical veclocity
%scatter (xr,yr)

t = 1.5e-14; %time interval that captures line
xd = vx*t; %displacement in x during one time interval
yd = vy*t; %displacement in y during one time interval

for p = 1:1:numOfStep
    xr = xr+xd;
    yr = yr+yd;
%%%%%%%%%Calculate average temperature of all particles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v = sqrt(vx.^2+vy.^2);
    TParticles = (mn*v.^2)/(C.kb);%Tempearture of individual particles
    Tave (p) = sum(TParticles)/numOfAtom;%Average temperature of all particles
    figure (1)
    plot (Tave)
    xlim ([0, numOfStep])
    ylim ([0, 500])
    xlabel("Number of steps (1.5e-14s/step)")
    ylabel("Temperature (K)")
    title("average temperature over time")
    
    %Define the left&right wrap-around
    xrp(xr >2e-7) = - (2e-7 - xrp(xr >2e-7));% changing previous point to prevent line drawn across canvas
    xr(xr > 2e-7) = xr(xr > 2e-7)-(2e-7);
    xrp(xr <0)    = 2e-7 - xrp(xr <0);% changing previous point to prevent line drawn across canvas
    xr(xr < 0)    = xr(xr < 0 )+(2e-7);
    
    %Define the specular top&bottom
    yd(yr > 1e-7) = - yd(yr > 1e-7);
    yr(yr > 1e-7) = (1e-7)-(yr(yr > 1e-7)-(1e-7));
    yd(yr < 0) = -yd(yr < 0 );
    yr(yr < 0) = -yr(yr < 0);
    figure (2)
    plot([xrp(1), xr(1)], [yrp(1), yr(1)], 'r')
    plot([xrp(2), xr(2)], [yrp(2), yr(2)], 'b')
    plot([xrp(3), xr(3)], [yrp(3), yr(3)], 'k')
    plot([xrp(4), xr(4)], [yrp(4), yr(4)], 'g')
    plot([xrp(5), xr(5)], [yrp(5), yr(5)], 'y')
    plot([xrp(6), xr(6)], [yrp(6), yr(6)], 'c')
    xlabel("Semiconductor Dimension (m)")
    ylabel("Semiconductor Dimension (m)")
    title ("Particles Trajectory")
    xlim ([0, 2e-7])
    ylim([0,1e-7])
    grid on
    hold on
    pause(0.05)
    
    xrp = xr;
    yrp = yr;
end
end

