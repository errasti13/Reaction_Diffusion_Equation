%-------------------------------------------------------------------------%
% Solving Adimensional Chemical Reaction-Diffusion Equation               %
%-------------------------------------------------------------------------%
% This MATLAB script solves the adimensional chemical reaction-diffusion  %
% equation for varying values of the Peclet number. The equation describes%
% the concentration profile of a chemical species undergoing reaction and %
% diffusion in a one-dimensional domain. The script uses the ode45 solver %
% to integrate the equation and produces plots of the concentration and   %
% normalized reaction rate.                                               %
%                                                                         %
% Author: Jon Errasti                                                     %
% Date: April 2022                                                        %
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

% Clear the workspace and command window
clear
clc

% Initialize parameters and constants
initialPeclet = 0.056;
targetError = 10^-3;
increment = 0.0001;
beta = 25;
Q = 1;
Da = 1;

% Loop for solving with increasing Peclet values
currentPeclet = initialPeclet;
count = 0;
while true
    currentPeclet = currentPeclet + increment;
    count = count + 1;

    % Define the differential equation and conversion to vector field
    syms y(t)
    diff_eq = diff(y, 2) == currentPeclet * diff(y, 1) + Da * y * exp(-beta * y * (1 + Q) / (1 + Q * (1 - y)));
    [V] = odeToVectorField(diff_eq);
    
    % Convert the vector field to a function handle
    M = matlabFunction(V, 'vars', {'t', 'Y'});
    
    % Set ODE solver options
    options = odeset('MaxStep', 0.01);
    
    % Solve the differential equation using ode45
    sol = ode45(M, [100, -100], [0.00001, 0], options);

    % Calculate the error and check if it meets the target threshold
    finalConcentration = sol.y(1, end);
    error = abs(finalConcentration(1) - 1);
    disp(['Iteration:', num2str(count), ' - Error:', num2str(error), ' - Peclet:', num2str(currentPeclet)]);
    
    if error <= targetError
        break;
    end
end

% Calculate normalized reaction rate and temperature adimensional
omega = sol.y(1, :) .* exp(-beta .* (sol.y(1, :) .* (1 + Q) ./ (1 + Q * (1 - sol.y(1, :)))));
omega = omega / max(omega);
T = (1 - sol.y(1, :)) * Q + 1;

%% Plotting
% Plot for normalized reaction rate
figure()
plot(sol.y(1, :), omega)
xlabel('$\mathcal{C}$', 'Interpreter', 'latex')
ylabel('$\frac{\Omega}{\Omega_{\max}}$', 'Interpreter', 'latex')
axis([0, 1, 0, 1])
legend('$\frac{\Omega}{\Omega_{\max}}$', 'Interpreter', 'latex', 'Location', 'SouthWest');

% Plot for concentration and temperature adimensional
figure()
plot(sol.x, sol.y(1, :))
ylabel('Concentración adimensional de combustible')
axis([-100, 100, 0, 1])
hold on
yyaxis right
plot(sol.x, T)
xlabel('$\xi$', 'Interpreter', 'latex')
ylabel('$\mathcal{T}$''$\frac{\Omega}{\Omega_{\max}}$', 'Interpreter', 'latex')
title(['Solucion para Pe=', num2str(currentPeclet)])
plot(sol.x, omega)
legend('$\mathcal{C}$', '$\mathcal{T}$', '$\frac{\Omega}{\Omega_{\max}}$', 'Interpreter', 'latex', 'FontSize', 18, 'Location', 'SouthWest');
