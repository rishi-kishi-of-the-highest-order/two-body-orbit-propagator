%-------------------------------------------------------------------%
% Rishi Kishore
%-------------------------------------------------------------------%

close all;
clear;
clc;

%-------------------------------------------------------------------%
% Problem 1
% Inputs:
%   position and velocity vector in ECI coords,
%   initial and later time,
%   relative and abs tolerance for iterative procedure
% Outputs:
%   input data with units
%   classical orbital elements
%   DCM converting from ECI to Perifocal Coords at time t1
%   position and velocity vectors in perifocal coords
%   True and Eccentric Anomalies at time t1
%   Time of Periapsis Passage
%   True and Eccentric Anomalies at time t2
%   Position and velocity vectors in Perifocal Coords at time t2
%   Position and velocity vectors in ECI coords at time t2
%   outputs must have labels and units
%-------------------------------------------------------------------%

%--------------- User-defined values for the problem ---------------%
r0_ECI = [5214.1, 6322.9, -3001.7]; % km (Position vector in ECI)
v0_ECI = [-4.1549, 3.1666, -4.5044]; % km/s (Velocity vector in ECI)
t1 = 0; % s (Initial time)
t2 = 12500; % s (Later time)
rel_tol = 1e-6; % Relative tolerance
abs_tol = 1e-6; % Absolute tolerance
mu = 398600; % km^3/s^2 (Gravitational parameter)

%-------------------------------------------------------------------%
% Part A: Displaying Name and Course Information
%-------------------------------------------------------------------%
fprintf("Rishi Kishore\n\n********** Orbital Mechanics Program **********\n");

%-------------------------------------------------------------------%
% Part B: Output Initial Condition Data
%-------------------------------------------------------------------%
fprintf("\nPart B: Output Initial Condition Data\n");
fprintf("Initial Position Vector: [%.5f, %.5f, %.5f] km\n", r0_ECI(1), r0_ECI(2), r0_ECI(3));
fprintf("Initial Velocity Vector: [%.5f, %.5f, %.5f] km/s\n", v0_ECI(1), v0_ECI(2), v0_ECI(3));
fprintf("\nInitial time: %.5f s\nLater time: %.5f s\n", t1, t2);
fprintf("\nRelative Tolerance: %.9f\nAbsolute Tolerance: %.9f\n", rel_tol, abs_tol);

%-------------------------------------------------------------------%
% Part C: Classical Orbital Elements (COE) Calculations
%-------------------------------------------------------------------%
fprintf("\nPart C: Calculating Classical Orbital Elements\n");

% Semi-major axis (a) using Vis-Viva Equation
r0_ECI_scalar = norm(r0_ECI); % km
v0_ECI_scalar = norm(v0_ECI); % km/s
a = -mu / (2 * ((v0_ECI_scalar^2) / 2 - mu / r0_ECI_scalar)); % km
fprintf("a = %.5f km\n", a);

% Eccentricity (e)
h_ECI = cross(r0_ECI, v0_ECI); % Specific angular momentum vector
h_ECI_scalar = norm(h_ECI);
e = sqrt(1 - ((h_ECI_scalar^2) / mu) / a); % Magnitude of eccentricity
fprintf("e (magnitude) = %.5f\n", e);
e_vector = (cross(v0_ECI, h_ECI) / mu) - (r0_ECI / r0_ECI_scalar); % Eccentricity vector
fprintf("e (vector) = [%.5f, %.5f, %.5f]\n", e_vector(1), e_vector(2), e_vector(3));

% Inclination (i)
i = acos(h_ECI(3) / h_ECI_scalar); % radians
fprintf("i = %.5f rad\n", i);

% Argument of Periapsis (ω)
K_north_pole = [0, 0, 1]; % North pole unit vector
n_hat = cross(K_north_pole, h_ECI) / norm(cross(K_north_pole, h_ECI));
if dot(e_vector, K_north_pole) > 0
    w = acos(dot(e_vector, n_hat) / e); % Northern Hemisphere
    fprintf("Periapsis in Northern Hemisphere: ω = %.5f rad\n", w);
else
    w = 2*pi - acos(dot(e_vector, n_hat) / e); % Southern Hemisphere
    fprintf("Periapsis in Southern Hemisphere: ω = %.5f rad\n", w);
end

% Right Ascension of the Ascending Node (RAAN, Ω)
I_hat = [1, 0, 0]; % First Point of Aries unit vector
J_hat = [0, 1, 0]; % Another unit vector for calculation
i_dot_n = dot(I_hat, n_hat);
j_dot_n = dot(J_hat, n_hat);
if i_dot_n > 0
    capOmega = atan(j_dot_n / i_dot_n);
else
    capOmega = pi + atan(j_dot_n / i_dot_n); % Adjustment for negative case
end
fprintf("Right Ascension of the Ascending Node (RAAN) = %.5f rad\n", capOmega);

%-------------------------------------------------------------------%
% Part D: DCM for ECI to Perifocal Coordinates
%-------------------------------------------------------------------%
fprintf("\nPart D: Calculating the DCM for ECI to Perifocal\n");

% Rotation Matrices
C1 = rotz(w); % Rotation about z-axis (ω)
C2 = rotx(i); % Rotation about x-axis (i)
C3 = rotz(capOmega); % Rotation about z-axis (RAAN)

% DCM that transforms ECI to Perifocal
C_EP = C1 * C2 * C3;
C_PE = C_EP'; % Inverse DCM (ECI to Perifocal)

fprintf("\nDCM to convert from ECI to Perifocal:\n");
fprintf("| %.6f   %.6f   %.6f |\n| %.6f   %.6f   %.6f |\n| %.6f   %.6f   %.6f |\n", ...
    C_EP(1, 1), C_EP(1, 2), C_EP(1, 3), C_EP(2, 1), C_EP(2, 2), C_EP(2, 3), C_EP(3, 1), C_EP(3, 2), C_EP(3, 3));

%-------------------------------------------------------------------%
% Part E: Position and Velocity at t1 in Perifocal Coordinates
%-------------------------------------------------------------------%
fprintf("\nPart E: Converting r0 and v0 from ECI to Perifocal\n");

r0_perifocal = C_EP * r0_ECI'; % km
v0_perifocal = C_EP * v0_ECI'; % km/s

fprintf("\nInitial Position Vector in Perifocal:\n [%.5f, %.5f, %.5f] km\n", r0_perifocal(1), r0_perifocal(2), r0_perifocal(3));
fprintf("Initial Velocity Vector in Perifocal:\n [%.5f, %.5f, %.5f] km/s\n", v0_perifocal(1), v0_perifocal(2), v0_perifocal(3));

%-------------------------------------------------------------------%
% Part F: True and Eccentric Anomalies at t1
%-------------------------------------------------------------------%
fprintf("\nPart F: True and Eccentric Anomalies at t1\n");

theta1 = acos((a * (1 - e^2)) / (r0_ECI_scalar * e) - 1 / e);
E1 = 2 * atan(tan(theta1 / 2) * sqrt((1 - e) / (1 + e))); % Eccentric Anomaly

fprintf("True Anomaly at initial time: %.5f rad\nEccentric Anomaly at initial time: %.5f rad\n", theta1, E1);

%-------------------------------------------------------------------%
% Part G: Time of Periapsis Passage
%-------------------------------------------------------------------%
fprintf("\nPart G: Time of Periapsis Passage (T_0)\n");

T0 = -(E1 - e * sin(E1)) * sqrt((a^3) / mu) + t1;
fprintf("T_0 = %.5f seconds\n", T0);

%-------------------------------------------------------------------%
% Part H: True and Eccentric Anomalies at t2
%-------------------------------------------------------------------%
fprintf("\nPart H: True and Eccentric Anomalies at t2\n");

% Newton-Raphson method for Eccentric Anomaly (E) at t2
E_new = sqrt(mu / a^3) * (t2 - T0); % Initial guess
E_old = 1e6 * E_new; % Large number for iteration

while abs(E_new - E_old) > abs_tol && abs((E_new - E_old) / E_old) > rel_tol
    E_old = E_new;
    E_new = E_old - (E_old - e * sin(E_old) - sqrt(mu / a^3) * (t2 - T0)) / (1 - e * cos(E_old));
end

theta2 = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E_new / 2)); % True Anomaly at t2

fprintf("True Anomaly at final time: %.5f rad\nEccentric Anomaly at final time: %.5f rad\n", theta2, E_new);

%-------------------------------------------------------------------%
% Rotation Matrices Functions
%-------------------------------------------------------------------%
function R = rotx(theta)
    % Rotation about x-axis
    R = [1, 0, 0;
         0, cos(theta), sin(theta);
         0, -sin(theta), cos(theta)];
end

function R = roty(theta)
    % Rotation about y-axis
    R = [cos(theta), 0, sin(theta);
         0, 1, 0;
         -sin(theta), 0, cos(theta)];
end

function R = rotz(theta)
    % Rotation about z-axis
    R = [cos(theta), sin(theta), 0;
        -sin(theta), cos(theta), 0;
         0, 0, 1];
end

