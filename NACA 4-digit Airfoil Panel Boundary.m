% SAMPLE CODE, THE START OF THE PANEL METHOD CODE, FOR AIRFOIL GEOMETRY

%==========================================================================
%
% VORTEX PANEL METHOD CODE FOR MAE 424: AERODYNAMICS
%
% 2nd-order panel method, i.e. uses panels with linearly-varying strength.
% First part of code creates (x,y) coordinates for NACA 4-digit airfoils.
% Second part calculates the lift coefficient, etc.
% An overall 'for' loop over k produces variations with angle of attack.
%
% Originally created for MATLAB under the supervision of Dr. Cyrus Madnia
% Revised in April 2015 and 2017 by Dr. Matthew Ringuette
% 
% Based on the panel method code published in:
% Kuethe & Chow "Foundations of Aerodynamics: Bases of
% Aerodynamic Design" 1998, John Wiley & Sons
% 
% That code based on: Stevens, W. A., Goradia, S. H., & Braden, J. A.,
% "Mathematical Model for Two-Dimensional Multi-Component Airfoils in
% Viscous Flow," NASA CR-1843, 1971.
%
% The airfoil geometry information is based on both Kuethe & Chow and:
% Ladson, C. L., Brooks, Jr., C. W., & Hill, A. S., "Computer Program to
% Obtain Ordinates for NACA Airfoils," NASA TM 4741, 1996.
%                                                               
%==========================================================================

clear all
clc

%==========================================================================
% Part I: user input (airfoil geometry, # panels)
%==========================================================================

pmxx = input('Enter desired NACA 4-digit designation, e.g. 0012: ');
disp(' ');
Np = input('Enter number of vortex panels (even integer): ');
disp(' ');

%==========================================================================
% Part II: creation of NACA 4-digit airfoil panel geometry
%
% Calculates the coordinates of the airfoil panel boundary points
% (xfoil,yfoil), and coordinates of the mean camber line points (xc,yc)
%==========================================================================

disp('Creating (x,y) coordinates of airfoil boundary points, camber line...');
disp(' ');

% Isolate key parameters of airfoil geometry
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pm    =	floor(pmxx/100);
xx    =	pmxx-pm*100; % Thickness in percent chord
p     = floor(pm/10);
m     = pm-p*10;
p     = p/100; % Rescale, max camber in chord lengths
m     = m/10; % Rescale, max camber location in chord lengths

% Compute the thickness distribution
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%

% Inititalize values for constants
a0 = 0.2969;
a1 = -0.1260;
a2 = -0.3516;
a3 = 0.2843;
a4 = -0.1015;

% Compute angle partition magnitude based on number of panels in input, Np
angle = (2*pi)/Np;

% Compute individual angles and xc and yt for each point
for i = 1:Np+1
    B(i) = angle*(i - 1);
    xc(i) = (0.5) + (0.5)*cos(B(i));
end

for j = 1:Np+1
    if j < Np/2+1
        yt(j) = -((xx)/20)*((a0*(xc(j).^(0.5))) + (a1*(xc(j))) + (a2*(xc(j).^2)) + (a3*(xc(j).^3)) + (a4*(xc(j).^4)));
    else
        yt(j) = ((xx)/20)*((a0*(xc(j).^(0.5))) + (a1*(xc(j))) + (a2*(xc(j).^2)) + (a3*(xc(j).^3)) + (a4*(xc(j).^4)));
    end
end

% Add the mean camber line
% ~~~~~~~~~~~~~~~~~~~~~~~~
%

% Camber line points computed
for n = 1:Np+1
    if xc(n) <= m
        yc = (p/(m^2))*((2*m*(xc)) - (xc).^2);
    else
        yc = (p/((1 - m)^2))*((1 - 2*m) + (2*m* (xc)) - (xc).^2);
    end
end

S = atan(yc/xc); % local inclination of the mean camber line

% New boundary points
xfoil = xc - yt.*sin(S); 
yfoil = yc + yt.*cos(S);

% Plot airfoil panel geometry to check the result
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
figure(1)
plot(xfoil(1:(Np/2+1)),yfoil(1:(Np/2+1)),'r'); % Lower airfoil boundary
hold on
plot(xfoil((Np/2+1):(Np+1)),yfoil((Np/2+1):(Np+1)),'g'); % Upper airfoil boundary
plot(xc,yc,'b'); % Chord line
daspect([1 1 1]);
xlabel('(x/c)');
ylabel('(y/c)');
title('Panel Boundary Points for any NACA 4-digit Airfoil');

