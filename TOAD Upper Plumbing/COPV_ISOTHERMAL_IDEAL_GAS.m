

close all;
%%%%%USING IDEAL GAS PV = mRT (assuming isothermal)%%%%%%%%%%%%%%%%%

%Pp = Pressure needed on propellant
%Vp = Volume prop tank, 
%Vo = Initial volume gas tank
%Pg = Gas tank pressue
%T0 = Initial temperature
%P0 = Initial pressure

function init_gas_mass = gasmass(Pp,Pg,Vp,T0, P0, pressurant)
if strcmp(pressurant, 'helium')
    R = 2079; % J/Mol*K
elseif strcmp(pressurant, 'nitrogen')
    R = 296.7;
end

init_gas_mass = ((Pp*Vp/(R*T0)) * (1/(1-(Pg/P0)))) * 10^6; %
end

sutton_example_1 = gasmass(3.40, 3.40, 0.180, 298, 14, 'nitrogen');
sutton_example_2 = gasmass(3.40, 3.40, 0.180, 298, 14, 'helium');

%tank_volume = input('Enter tank volume here: ') * (2.54 * 1e-2)^3;
tank_volume = 200 * (2.54 * 1e-2)^3;
%init_tank_pressure = input('Enter init tank pressure here: ');
init_tank_pressure = 4500 * 6894.76;
%throat_d = input('Enter throat diameter here: ') * (2.54 * 1e-2)^2;
%Cd = input("Enter discharge coefficent here: ");
Cd = .7;
seconds = input("Enter seconds here: ");
intervals = input("Enter intervals here: ");
%init_tank_density = (gasmass(3.40, 3.40, 0.180, 298, 14, 'helium')/tank_volume)


M = 4.003; %molecular weight helium
T0 = 293; %kelvin
R = 2079;
gamma = 1.67;
%Astar = (pi/4)*throat_d^2;
Astar = 0.05 * (2.54 * 1e-2)^2; %imperial
c0 = (gamma * (R*T0))^(1/2); 
init_tank_density = init_tank_pressure/(R * T0);

function tau_constant = time_constant(tank_volume, Cd, Astar, c0, gamma)
    tau_constant = (tank_volume/(Cd*Astar*c0))*((gamma+1)/2)^((gamma+1)/2*(gamma-1));
end

function tank_density = copv_density_calc(init_density, t, tau)
    tank_density = init_density * (exp(-(t/tau)));
end

function tank_pressure = copv_pressure_calc(init_pressure, t, tau)
    tank_pressure = init_pressure * (exp(-(t/tau)));
end

% function mass_flow_rate = throat_mdot(Cd, Astar, tank_pressure, tank_temp, gamma, R, pressure_down)
% 
%     mass_flow_rate = Cd *Astar * pressure_down *sqrt((2*gamma / (R * tank_temp * (gamma - 1)) * ((tank_pressure/pressure_down)^((gamma-1)/gamma)))*(1+((pressure_down/tank_pressure)^(gamma-1/gamma))));
% end
function mass_flow_rate = throat_mdot(Cd, Astar, tank_density, tank_temp, gamma, R, M)
    mass_flow_rate = Cd * Astar * tank_density * sqrt((gamma * R * tank_temp)) * ((2 / (gamma + 1))^((gamma + 1) / (2 * (gamma - 1))));
end

tau_1 = time_constant(tank_volume, Cd, Astar, c0, gamma)

% Time intervals
time_intervals = linspace(0, seconds, intervals);
blowdown_density_values = zeros(length(time_intervals));
% Initialize with the initial value
blowdown_pressure_values = zeros(length(time_intervals));
% Initialize with the initial pressure
blowdown_mdot = zeros(length(time_intervals));  % Initial mass flow rate

% Index for the loop
% Loop through the time intervals and calculate the new values
%pressure_down = 101325;
for index = 1:length(time_intervals)
    % Get the current time step
    current_time = time_intervals(index);
    
    % Calculate new density, pressure, and mass flow rate
    new_density = copv_density_calc(init_tank_density, current_time, tau_1);
    %fprintf("current_time: %f\n", current_time);

    % Correct the line here:
    new_pressure = copv_pressure_calc(init_tank_pressure, current_time, tau_1);
    %fprintf("new_pressure: %0.4f \n", new_pressure);

    new_mdot = throat_mdot(Cd, Astar, new_density, T0, gamma, R, M);

    % Append the new values to the arrays
    blowdown_density_values(index) = new_density;
    blowdown_pressure_values(index) = new_pressure;
    blowdown_mdot(index) = new_mdot;
end

% Plot the results
figure;
plot(time_intervals, blowdown_pressure_values);
title('Blowdown Pressure');
xlabel('Time (seconds)');
ylabel('Pressure (Pa)');

figure;
plot(time_intervals, blowdown_mdot);
title('mdot');
xlabel('Time (seconds)');
ylabel('mdot');

figure;
plot(time_intervals, blowdown_density_values);
title('Blowdown Density');
disp(blowdown_mdot(1));

%Display the final pressure values and time intervals
%disp(blowdown_pressure_values);
%blowdown_pressure_values(1)
% blowdown_pressure_values(15)
% blowdown_pressure_values(500)
% %disp(time_intervals);
%copv_pressure_calc(3.4, 0.375, 5.8374);
%disp(Astar);



%%%%NOTES%%%%%
%Add downstream press
%Fix mdot

