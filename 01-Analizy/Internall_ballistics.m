clear; close all;
%% Input data

%Performance required
thrust = 400;
vectoring = 0.25;
operation_time = 120;
ambient_pressure = 0.1;
firing_number = 15;

%Termodynamic properties interpolants of combustion products
load("Trade_off_mp\output_n2o\Propan.mat");

interp_of_ratio = x.output.oxfl' .* ones(1,size(x.output.oxfl,2));
interp_pressure = x.output.pressure / 10;
%[interp_of_ratio, interp_pressure] = meshgrid(interp_of_ratio, interp_pressure);

temperature_matrix = x.output.temperature;
kappa_matrix = x.output.gamma;
molar_mass_matrix = x.output.mw;

%set interpolants
Combustion_chamber.temperature_interp = griddedInterpolant({x.output.oxfl', interp_pressure(1,:)},...
    temperature_matrix, 'linear', 'nearest');
Combustion_chamber.kappa_interp = griddedInterpolant({x.output.oxfl', interp_pressure(1,:)},...
    kappa_matrix, 'linear', 'nearest');
Combustion_chamber.molar_mass_interp = griddedInterpolant({x.output.oxfl', interp_pressure(1,:)},...
    molar_mass_matrix, 'linear', 'nearest');

%Combustion chamber conditions assumed
Combustion_chamber.pressure = 1.3;
Combustion_chamber.of = 9;

%% Calculations

%Actual thermodynamic parameters in chamber
Combustion_chamber.temperature = Combustion_chamber.temperature_interp(Combustion_chamber.of,...
    Combustion_chamber.pressure);
Combustion_chamber.kappa = Combustion_chamber.kappa_interp(Combustion_chamber.of,...
    Combustion_chamber.pressure);
Combustion_chamber.molar_mass = Combustion_chamber.molar_mass_interp(Combustion_chamber.of,...
    Combustion_chamber.pressure);
Combustion_chamber.gas_constant = 8314 / Combustion_chamber.molar_mass;

Combustion_chamber.thrust_coef = sqrt(2 * Combustion_chamber.kappa ^ 2 / (Combustion_chamber.kappa...
    - 1) * (2 / (Combustion_chamber.kappa + 1)) ^ ((Combustion_chamber.kappa + 1) /...
    (Combustion_chamber.kappa - 1)) * (1 - (ambient_pressure / Combustion_chamber.pressure)));
Combustion_chamber.nozzle_throat_diameter = sqrt(4 * thrust / pi / Combustion_chamber.pressure /...
    Combustion_chamber.thrust_coef);
Combustion_chamber.nozzle_throat_area = pi * Combustion_chamber.nozzle_throat_diameter ^ 2 / 4;
Combustion_chamber.c_star = sqrt(Combustion_chamber.kappa * Combustion_chamber.gas_constant *...
    Combustion_chamber.temperature) / (Combustion_chamber.kappa * sqrt( (2 /...
    (Combustion_chamber.kappa+1) ) ^ ((Combustion_chamber.kappa + 1) / (Combustion_chamber.kappa - 1) )));
Combustion_chamber.nozzle_effective_exhaust_velocity = Combustion_chamber.c_star * Combustion_chamber.thrust_coef;
Combustion_chamber.sonic_velocity = sqrt(Combustion_chamber.kappa * Combustion_chamber.gas_constant *...
    Combustion_chamber.temperature);
Combustion_chamber.nozzle_exhaust_mach = Combustion_chamber.nozzle_effective_exhaust_velocity /...
    Combustion_chamber.sonic_velocity;
Combustion_chamber.nozzle_area_ratio = 1 / Combustion_chamber.nozzle_exhaust_mach *...
    ((1 + (Combustion_chamber.kappa - 1) * Combustion_chamber.nozzle_exhaust_mach ^ 2 / 2) / (1 +...
    (Combustion_chamber.kappa - 1) / 2)) ^ ((Combustion_chamber.kappa + 1) / 2 / (Combustion_chamber.kappa - 1));
Combustion_chamber.nozzle_exit_area = Combustion_chamber.nozzle_throat_area * Combustion_chamber.nozzle_area_ratio;
Combustion_chamber.nozzle_exit_diameter = sqrt(4 * Combustion_chamber.nozzle_exit_area / pi);








