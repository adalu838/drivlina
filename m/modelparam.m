% modelparam.m

%

% Defines workspace variables which

% are used by the simulink model.

%

% Aims to resemble a passenger car

%

%Konstanter som ej var definerade

i_t = 3.4;
J_t = 0.13;
m = 1500;
k = 1000;
c = 1;
c_r1 = 0.147;



% Environment

% -----------

% road topography

position = [0 1e6]; slope = [0 0]; altitude = [0 0];

% initial speed [km/h]

init_v = 30;

% Input signals

lowlevel = -0.02;

highlevel = 0.2;

t1 = 12; t2 = 17;

InputT = [0 t1-eps t1 t2-eps t2]';

InputU = [lowlevel lowlevel highlevel highlevel lowlevel]';

InputG = [1 1 1 1 1]';



%

% Chassis

% -------

% mass [kg]

mass = 1.5e3;

% gravity constant [m/s^2]

g = 9.81;

% wheel inertia [kg m^2]

J_w = 0.6;

% wheel radius [m]

r_w = 0.30;

% air drag coefficient [-]

c_w = 0.30;

% maximum vehicle cross section area [m^2]

A_a = 1.5;

% Air density [kg/m^3]

Rho_a = 1.29;

% rolling resistance coefficient [-]

c_r = 0.015;

% initial drive shaft speed [rad/s]

init_d = init_v/r_w/3.6;

% brake gain [Nm]

% brake torque T is calculated as T = - brake_gain*B

% where B is the control signal [0,1]

brake_gain = 1.0e3;



%

% Drive shafts

% ------------

% stiffness [Nm/rad]

k_d = 1e3; %1e3;

% damping [Nm s/rad]

c_d = 1;

% final drive initial speed [rad/s]

init_f = init_d;



%

% Final drive

% -----------

% inertia [kg m^2]

J_f = 0.1;

% ratio [-]

i_f = 3.4;

% Damping [Nm s/rad]

b_f = 0.1;



%

% Propeller shaft

% ---------------

% stiffness [Nm/rad]

k_p = 0.5e4; %5e3;

% damping [Nm s/rad]

c_p = 15;

% gearbox initial speed [rad/s]

init_t = init_f*i_f;



%

% Gearbox

% -------

% Damping [Nm s/rad]

b_t = 0.1;

% gear number to ratio and inertia tables

gearNumber = [1 2 3 4 5];

gearRatio = [3.4 2.2 1.5 1.1 0.9];

gearInertia = [0.13 0.09 0.07 0.05 0.04];



%

% Clutch

% ------

% stiffness [Nm/rad]

k_c = 10e3; %10e3;

% damping [Nm s/rad]

c_c = 30;

% initial engine speed [rad/s]

init_e = init_t*gearRatio(2);



%

% ICE

% ----

% engine inertia [kg m^2]

J_e = 0.2;

% time delay [s]

tau_d = 1e-3;

% time constant [s]

tau_e = 0.1;

% max torque [Nm]

T_max = 400;

% max speed [rad/s]

N_max = 7e3*pi/30;

% idle speed [rad/s]

N_idle = 8.5e2*pi/30;

% idle controller proportional gain

IdleControllerP = 0.8;

% idle controller integral gain

IdleControllerI = 0.2;

% initial torque

init_T = 0;%0.5*c_w*A_a*Rho_a*r_w^3*init_d^2 + mass*g*r_w*c_r;



%

% Longitudinal force

% ------------------

% Longitudinal force, pure slip, Pacjeka (2002), p. 187

% all lambda = 0, all \zeta = 1, all \eps = 0

% dfz = 0; pHx1 = 0; pHx2 = 0; pEx4 = 0;

Fz = 0.5*mass*g;

Kxk = 12*Fz;

Cx = 1.2;

muX = 0.8;

Dx = muX*Fz;

Bx = Kxk/(Cx*Dx);

Ex = 0.8;

%Fxo = Dx*sin(Cx*atan(Bx*kappaX-Ex.*(Bx*kappaX-atan(Bx*kappaX))))+Svx;



% NOTE

% Pacejka uses the SAE definiton of long. slip:

%  -\frac{\Omega_0-\Omega}{\Omega_0}

% The definition used in TSFS05 is

%  -\frac{\Omega_0-\Omega}{\Omega}

% where \Omega_0 = V_x/r_e