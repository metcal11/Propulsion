


mdot_IPA = input("Input IPA Mdot (lb/s): "); % Get mass flow rate
mdot_LOX = input("Input LOX Mdot (lb/s): "); % Get mass flow rate
multiplier = input("Multiplier: ");

mdot_LOX = mdot_LOX * multiplier;
mdot_IPA = mdot_IPA * multiplier;

%thickness parameters
y_strength = 40000; %PSI for 6061 T-6
internal_pressure = input("MEOP tank pressure (psi): ");
input_mos = input("Insert MOS: ");


%%ellipsoid = 4/3pi a*a*b (a is radius, b is height. ratio should be
%%sqrt(2)

flight_time = input("Flight time (s): "); % Get flight duration
%of_ratio = input("Input O/F ratio: ");
%x = of_ratio + 1;

%mdot_IPA = mdot / x
%mdot_LOX = 1.2 * mdot_IPA

density_ipa = 6.5; % lbs per gallon
density_lox = 9.52;

ipa_mass = mdot_IPA * flight_time; % Total ipa mass
ipa_vol = ipa_mass / density_ipa; % Volume in gallons
ipa_vol_cbft = ipa_vol / 7.48052; % Convert to cubic feet

lox_mass = mdot_LOX * flight_time; % Total lox mass
lox_vol = lox_mass / density_ipa; % Volume in gallons
lox_vol_cbft = lox_vol / 7.48052; % Convert to cubic feet

fprintf("IPA Tank Volume: %.4f cubic feet\n", ipa_vol_cbft);
fprintf("IPA mass: %.02f lb\n", ipa_mass);

fprintf("LOX Tank Volume: %.4f cubic feet\n", lox_vol_cbft);
fprintf("LOX mass: %.02f lb\n", lox_mass);

%%%%CHANGE THESE VARIABLES FOR TANK SIZING%%%%%%%%%%%
max_radius_ipa = 2/3; % ft
min_radius_ipa = 0.25; % ft
max_height_ipa = 5; % ft

max_radius_lox = 0.66; % ft
min_radius_lox = 0.25; % ft
max_height_lox = 5; % ft

% Define Radius Range
ipa_radius = linspace(min_radius_ipa, max_radius_ipa, 100);
lox_radius = linspace(min_radius_lox, max_radius_lox, 100);
% Compute Heights

height_ipa = ipa_vol_cbft ./ (pi * ipa_radius.^2);
height_lox = lox_vol_cbft ./ (pi * lox_radius.^2);

cap_height_ipa = ipa_radius/sqrt(2);
cap_height_lox = lox_radius/sqrt(2);
new_cyl_ipa_height = (ipa_vol_cbft) ./ ((pi * ipa_radius.^2)) - (4/3) * cap_height_ipa;
new_cyl_lox_height = (lox_vol_cbft) ./ ((pi * lox_radius.^2)) - (4/3) * cap_height_lox;
new_height_ipa = new_cyl_ipa_height + (2 * cap_height_ipa);
new_height_lox = new_cyl_lox_height + (2 * cap_height_lox);

%thickness
ipa_radius_in = ipa_radius * 12;
desired_thickness_thin = ((internal_pressure * ipa_radius_in) * (input_mos + 1) / (y_strength)); %took 1/2 out, Grant told me to multiply by 2 but I dont believe him.
desired_thickness_thick = (sqrt((internal_pressure * (ipa_radius_in.^2) + ((ipa_radius_in.^2) * y_strength)) / (y_strength - (internal_pressure * input_mos) - internal_pressure)) - ipa_radius_in); %took 1/2 out, %took 1/2 out, Grant told me to multiply by 2 but I dont believe him.

% Ensure heights don’t exceed max_height
valid_indices_ipa = new_height_ipa <= max_height_ipa; % Logical mask for valid region
radius_valid_ipa = ipa_radius(valid_indices_ipa); % Keep only valid radii
height_valid_ipa = new_height_ipa(valid_indices_ipa); % Keep only valid heights

valid_indices_lox = new_height_lox <= max_height_lox; % Logical mask for valid region
radius_valid_lox = lox_radius(valid_indices_lox); % Keep only valid radii
height_valid_lox = new_height_lox(valid_indices_lox); % Keep only valid heights

%%%%%%%%%%%%%%%%%%%%%%%%IPA Plot%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(ipa_radius, new_height_ipa, 'b', 'LineWidth', 2);
hold on;

% Constraint Lines
xline(max_radius_ipa, 'r--', 'Max Radius', 'LabelHorizontalAlignment', 'left');
yline(max_height_ipa, 'r--', 'Max Height', 'LabelVerticalAlignment', 'middle');

% Labels and Formatting
xlabel('Radius (ft)');
ylabel('Height (ft)');
title('Valid IPA Dimensions');
grid on;
legend('h = Vol / (\pi r^2)', 'Height Constraint', 'Location', 'NorthEast');
xlim([min_radius_ipa - 0.05, max_radius_ipa + 0.05]);
ylim([0 max_height_ipa + 2]);

% Fill the region between the curve and max_height
fill([radius_valid_ipa, fliplr(radius_valid_ipa)], [height_valid_ipa, ones(size(radius_valid_ipa)) * max_height_ipa], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

hold off;

%%%%%%%%%%%%%%%%%%%%%%%LOX Plot%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
plot(lox_radius, new_height_lox, 'b', 'LineWidth', 2);
hold on;

% Constraint Lines
xline(max_radius_lox, 'r--', 'Max Radius', 'LabelHorizontalAlignment', 'left');
yline(max_height_lox, 'r--', 'Max Height', 'LabelVerticalAlignment', 'middle');

% Labels and Formatting
xlabel('Radius (ft)');
ylabel('Height (ft)');
title('Valid LOX Dimensions');
grid on;
legend('h = Vol / (\pi r^2)', 'Height Constraint', 'Location', 'NorthEast');
xlim([min_radius_lox - 0.05, max_radius_lox + 0.05]);
ylim([0 max_height_lox + 2]);

% Fill the region between the curve and max_height
fill([radius_valid_lox, fliplr(radius_valid_lox)], [height_valid_lox, ones(size(radius_valid_lox)) * max_height_lox], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

hold off;


%%%%%%%%%%%%%%%%%%Tank Thickness%%%%%%%%%%%%%%%%%%%%%%%%%%



figure;
plot(ipa_radius_in, desired_thickness_thin, 'b', 'LineWidth', 2);
hold on;

% Define fill region
x_fill = [ipa_radius_in, fliplr(ipa_radius_in)]; % X values (mirrored)
y_fill = [desired_thickness_thin, ones(size(desired_thickness_thin)) * max(desired_thickness_thin) * 1.2]; % Y values (mirrored above the curve)

fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % blue fill above

xlabel("Tank Radius (in)");
ylabel("Min allowable thickness (in)");
title("Tank Thickness with Thin Walled Stress");
grid on;
hold off;

figure;
plot(ipa_radius_in, desired_thickness_thick, 'b', 'LineWidth', 2);
hold on;

% Define fill region
x_fill = [ipa_radius_in, fliplr(ipa_radius_in)]; % X values (mirrored)
y_fill = [desired_thickness_thick, ones(size(desired_thickness_thin)) * max(desired_thickness_thin) * 1.2]; % Y values (mirrored above the curve)

fill(x_fill, y_fill, 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none'); % blue fill above

xlabel("Tank Radius (in)");
ylabel("Min allowable thickness (in)");
title("Tank Thickness with Thick Walled Stress");
grid on;
hold off;





