% Orbital simulation for autonomous navigation in space project
% Two satellites in LEO — Sat A is the chaser (12U cubesat, GLOBAL-12)
% and Sat B is the target debris (Electron kick stage rocket body)

clear; clc; close all;

%% Constants
R = 6.371e6;        % mean Earth radius [m]
Re = 6378136.3;      % equatorial radius — used in J2 [m]
M_earth = 5.972e24;       % Earth mass [kg]
G = 6.67430e-11;    % gravitational constant [m^3/kg/s^2]
mu = G * M_earth;    % precalculated so I don't repeat G*M everywhere
J2 = 1.08262668e-3;  % J2 oblateness coefficient (WGS-84)
AU  = 1.496e11;       % 1 AU in metres, needed for SRP
P_sun = 4.56e-6;        % solar radiation pressure at 1 AU [N/m^2]

%% Satellite A properties
% Modelled after GLOBAL-12 payload (NORAD 49772)
% Box satellite: 0.5 x 0.5 x 1.0 m, 55 kg
satA.m = 55;
satA.Cd = 2.2;
satA.A_drag = 0.03;       % mean projected drag area: 0.5 x 0.5 face [m^2]
satA.Cr = 1.3;
satA.A_srp = 0.25;       % consistent with drag area [m^2]

% Inertia tensor — uniform solid box 
bx = 0.50; by = 0.50; bz = 1.00;
Ixx_A = (satA.m/12) * (by^2 + bz^2);   
Iyy_A = (satA.m/12) * (bx^2 + bz^2);   
Izz_A = (satA.m/12) * (bx^2 + by^2);   
satA.I = diag([Ixx_A, Iyy_A, Izz_A]);

% Residual magnetic dipole 
% Using ~0.04 A.m^2 total, consistent with a clean small satellite platform
satA.mag_dip = [0.03; 0.02; -0.01];   % [A.m^2] body frame

% Face geometry for SRP and aerodynamic torque calculations
satA.faces = {[1;0;0],[-1;0;0],[0;1;0],[0;-1;0],[0;0;1],[0;0;-1]};
satA.areas = [0.50, 0.50, 0.50, 0.50, 0.25, 0.25];   % [m^2]
satA.reflectivity = [1.35, 1.35, 1.35, 1.35, 1.50, 1.50];   % sides / top-bottom

% SRP centre of pressure offsets from CoM [m]
% Larger body = slightly larger possible CoP offset from panel asymmetry
% Using ~1-2% of body dimension as CoP offset estimate
satA.cp_srp = {[ 0.010;  0.008; -0.005], [-0.010; -0.008;  0.005], ...
                 [ 0.008;  0.010;  0.006], [-0.008; -0.010; -0.006], ...
                 [ 0.006; -0.005;  0.018], [-0.006;  0.005; -0.018]};

% Aerodynamic CoP offsets
satA.cp_aero = {[ 0.008;  0.006; -0.004], [-0.008; -0.006;  0.004], ...
                 [ 0.006;  0.008;  0.005], [-0.006; -0.008; -0.005], ...
                 [ 0.005; -0.004;  0.015], [-0.005;  0.004; -0.015]};


%% Satellite B properties
% Modelled after Electron Kick Stage R/B (NORAD 64230)
% Approximated as a cylinder: radius 0.475 m, length 1.2 m, 45 kg
satB.m = 45;
satB.Cd = 2.2;
satB.A_drag = 0.71;       % lateral projected area — tumbling body uses the large side [m^2]
satB.Cr = 1.15;        
satB.A_srp = 0.71;       

% Inertia — solid cylinder formulas
% Ixy = m*(3r^2 + L^2)/12, Iz = m*r^2/2
r_cyl = 0.475; L_cyl = 1.2;
Ixy_B = (satB.m/12) * (3*r_cyl^2 + L_cyl^2);
Iz_B  = (satB.m/2)  * r_cyl^2;
satB.I = diag([Ixy_B, Ixy_B, Iz_B]);

% Residual magnetic dipole for the rocket body
satB.mag_diop = [0.010; 0.008; -0.012];   % 0.0175[A.m^2]

% Face geometry — same 6-face box approximation
satB.faces        = {[1;0;0],[-1;0;0],[0;1;0],[0;-1;0],[0;0;1],[0;0;-1]};
satB.areas        = [0.71, 0.71, 0.71, 0.71, 0.71, 0.36];   % [m^2]
satB.reflectivity = [1.35, 1.35, 1.35, 1.35, 1.35, 1.35];
satB.cp_srp  = {[ 0.30;  0.12; -0.02], [-0.30; -0.04;  0.02], ...
                 [ 0.12;  0.30;  0.03], [-0.04; -0.30; -0.03], ...
                 [ 0.05; -0.03;  0.50], [-0.05;  0.02; -0.30]};
satB.cp_aero = {[ 0.28;  0.10; -0.02], [-0.28; -0.03;  0.02], ...
                 [ 0.10;  0.28;  0.03], [-0.03; -0.28; -0.03], ...
                 [ 0.04; -0.02;  0.45], [-0.04;  0.02; -0.28]};

%% Orbital initial conditions
% Both satellites start at the same reference position on the x-axis
% at the 408 km reference altitude — each then gets its own perigee velocity
Altitude_ref = 408e3;
r0 = [R + Altitude_ref; 0; 0];

% Orbit A — chaser (GLOBAL-12, NORAD 49772, inclination = 53 deg)
incA = deg2rad(53);
ra_A_m = R + 425e3;   rp_A_m = R + 420e3;
eccenA = (ra_A_m - rp_A_m) / (ra_A_m + rp_A_m);
semi_majorA = (ra_A_m + rp_A_m) / 2;
vel_perigeeA = sqrt(mu * (1+eccenA) / (semi_majorA * (1-eccenA)));
Rx_A = [1 0 0; 0 cos(incA) -sin(incA); 0 sin(incA) cos(incA)];
stateA0 = [r0; Rx_A * [0; vel_perigeeA; 0]];

% Orbit B — target (Electron Kick Stage, NORAD 64230, i = 59 deg)
incB = deg2rad(59);
ra_B_m = R + 414e3;   rp_B_m = R + 408e3;
eccenB = (ra_B_m - rp_B_m) / (ra_B_m + rp_B_m);
semi_majorB = (ra_B_m + rp_B_m) / 2;
v_perigeeB = sqrt(mu * (1+eccenB) / (semi_majorB * (1-eccenB)));
Rx_B  = [1 0 0; 0 cos(incB) -sin(incB); 0 sin(incB) cos(incB)];
stateB0 = [r0; Rx_B * [0; v_perigeeB; 0]];

%% Initial attitude for Satellite B
% Starting aligned with inertial frame (identity quaternion)
% with a small initial tumble rate — representative of post-separation spin
% First attempt used Euler angles — hit gimbal lock during the close approach - switched t Quaternions      
q_init = [1; 0; 0; 0];
w_init = deg2rad([0.5; 0.3; 0.5]);   % [rad/s]
stateAtt0 = [q_init; w_init];

%% Simulation time
Ta = 2*pi * sqrt(semi_majorA^3 / mu);
Tb = 2*pi * sqrt(semi_majorB^3 / mu);
N = 1;   % Number of orbits, set to 1, change to >1, e.g. 10 to see orbital energy decay (value of 10 used in Dissertation)
tspan = [0, N * Tb];
fprintf('Orbital periods:\n  A: %.2f min\n  B: %.2f min\n', Ta/60, Tb/60);


%% Orbital integration
% Tight tolerances to keep position error small over one orbit
opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
[tA, sA] = ode45(@(t,s) OrbitsGeneral(t, s, satA, J2, Re, mu), tspan, stateA0, opts);
[tB, sB] = ode45(@(t,s) OrbitsGeneral(t, s, satB, J2, Re, mu), tspan, stateB0, opts);
xA = sA(:,1); yA = sA(:,2); zA = sA(:,3);
xB = sB(:,1); yB = sB(:,2); zB = sB(:,3);

%% interpolants for the attitude ODE
orb_A = @(t) interp1(tA, sA, t, 'linear')';
orb_B = @(t) interp1(tB, sB, t, 'linear')';
%% Sat B attitude integration
% No control — just disturbance torques driving the tumble
opts_att = odeset('RelTol', 1e-8, 'AbsTol', 1e-11);
[tAtt_B, sAtt_B] = ode45(@(t,s) AttitudeDynamics(t, s, orb_B, satB), ...
                          tspan, stateAtt0, opts_att);
q_hist_B = sAtt_B(:, 1:4);
w_hist_B = sAtt_B(:, 5:7);

%% Sat A pointing attitude (computed geometrically, not integrated)
% Sat A always points its x-body axis toward Sat B
% Interpolate Sat A state onto Sat B time grid first because ode45 gives
% different numbers of steps for each integration
sA_on_tB = interp1(tA, sA, tB, 'linear');
q_hist_A = zeros(length(tB), 4);
w_hist_A = zeros(length(tB), 3);
for i = 1:length(tB)
    r_A = sA_on_tB(i, 1:3)';
    r_B = [xB(i); yB(i); zB(i)];
    v_A = sA_on_tB(i, 4:6)';
    v_B = [sB(i,4); sB(i,5); sB(i,6)];
    dv = r_B - r_A;
    dv_n = norm(dv);
    if dv_n < 1e-3
        q_hist_A(i,:) = [1 0 0 0]; continue   % satellites coincide — skip
    end
    x_body = dv / dv_n;   % body x-axis points toward target
    z_tmp = cross(r_A, x_body);
    if norm(z_tmp) > 1e-6
        z_body = z_tmp / norm(z_tmp);   % orbit normal
    else
        z_body = [0; 0; 1];    % fallback if vectors are parallel
    end
    y_body = cross(z_body, x_body);     % complete right-hand frame
    C_IB = [x_body, y_body, z_body];
    q_hist_A(i,:) = dcm2quat(C_IB);
    if i > 1
        % Angular velocity from transport theorem: w = r_hat x r_hat_dot
        r_rel = r_B - r_A;   v_rel = v_B - v_A;
        rn = norm(r_rel);
        xdot = (v_rel - (dot(r_rel,v_rel)/rn^2)*r_rel) / rn;
        w_hist_A(i,:) = (C_IB' * cross(x_body, xdot))';
    end
end

%% 3D Animation
figure('Color', 'white', 'Position', [100, 100, 1200, 900]);
[Xs,Ys,Zs] = sphere(100);
earthSurf = surf(Xs*R, Ys*R, Zs*R, 'FaceColor', [0.2 0.3 1], ...
    'EdgeColor', 'white', 'FaceAlpha', 0.7);
hold on;
orbitA_line = plot3(xA, yA, zA, 'r', 'LineWidth', 1.5);
orbitB_line = plot3(xB, yB, zB, 'g', 'LineWidth', 1.5);
axis equal; grid on;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('Orbital simulation — J2, drag, SRP, magnetic and aero torques');
view(45, 30);
markerA = plot3(xA(1),yA(1),zA(1),'ro','MarkerFaceColor','r','MarkerSize',8);
markerB = plot3(xB(1),yB(1),zB(1),'mo','MarkerFaceColor','m','MarkerSize',8);
sc = 2e6;                                                                       % scale factor so attitude arrows are visible at orbital scale on the 3d model
CbiA0 = quat2dcm(q_hist_A(1,:)');
arrowA = quiver3(xA(1),yA(1),zA(1), sc*CbiA0(1,1), sc*CbiA0(2,1), sc*CbiA0(3,1), ...
    'c', 'LineWidth', 2, 'MaxHeadSize', 0.3, 'AutoScale', 'off');
CbiB0 = quat2dcm(q_hist_B(1,:)');
arrowB = quiver3(xB(1),yB(1),zB(1), sc*CbiB0(1,1), sc*CbiB0(2,1), sc*CbiB0(3,1), ...
    'y', 'LineWidth', 2, 'MaxHeadSize', 0.3, 'AutoScale', 'off');
% Legend placed after all handles are defined so nothing gets auto-named
legend([earthSurf, orbitA_line, orbitB_line, markerA, markerB, arrowA, arrowB], ...
       {'Earth', 'Global-12', 'Electron Kick Stage', ...
        'Global-12', 'Electron Kick Stage', ...
        'Global-12 Orientation', 'Electron Kick Stage Orientation'}, ...
       'Location', 'northeast', 'TextColor', 'white');
n = length(tB); speedFactor = 500;   % run 500x faster than real time 
% change to higher for faster results, otherise waiting 10-20 seconds just
% for model to play before graphs show up (or just cut this section and
% repaste it once graphs obtained)
for i = 1:n
    set(markerA, 'XData', xA(i), 'YData', yA(i), 'ZData', zA(i));
    set(markerB, 'XData', xB(i), 'YData', yB(i), 'ZData', zB(i));
    CbiA = quat2dcm(q_hist_A(i,:)');
    set(arrowA, 'XData', xA(i), 'YData', yA(i), 'ZData', zA(i), 'UData', sc*CbiA(1,1), 'VData', sc*CbiA(2,1), 'WData', sc*CbiA(3,1));
    CbiB = quat2dcm(q_hist_B(i,:)');
    set(arrowB, 'XData', xB(i), 'YData', yB(i), 'ZData', zB(i),'UData', sc*CbiB(1,1), 'VData', sc*CbiB(2,1), 'WData', sc*CbiB(3,1));
    if mod(i,5)==0, drawnow; % only refresh every 5 steps (otherwise becomes buggy)
    end   
    if i < n, pause((tB(i+1)-tB(i))/speedFactor); end
end
fprintf('Animation complete!\n');

%% Torque histories
% Re-evaluating magnetic, SRP and aero torques along each satellite's trajectory
% so can be plotted and observed for verification 
nA = length(tB);
nB = length(tAtt_B);
tau_mag_A  = zeros(nA,3); tau_srp_A  = zeros(nA,3); tau_aero_A = zeros(nA,3);
for i = 1:nA
    r_s = [xA(i); yA(i); zA(i)];
    v_s = sA_on_tB(i, 4:6)';
    q = q_hist_A(i,:)';
    C = quat2dcm(q);
    B = earthMagneticField(r_s);
    tau_mag_A(i,:) = cross(satA.mag_dip, C'*B)';
    tau_srp_A(i,:) = srpTorque(r_s, q, getSunPosition(tB(i)), satA)';
    tau_aero_A(i,:) = aeroDragTorque(r_s, q, v_s, satA)';
end
tau_mag_B  = zeros(nB,3); tau_srp_B  = zeros(nB,3); tau_aero_B = zeros(nB,3);
for i = 1:nB
    sv  = orb_B(tAtt_B(i));
    r_s = sv(1:3);   v_s = sv(4:6);
    q = sAtt_B(i, 1:4)';
    C  = quat2dcm(q);
    B = earthMagneticField(r_s);
    tau_mag_B(i,:) = cross(satB.mag_diop, C'*B)';
    tau_srp_B(i,:) = srpTorque(r_s, q, getSunPosition(tAtt_B(i)), satB)';
    tau_aero_B(i,:) = aeroDragTorque(r_s, q, v_s, satB)';
end
tau_total_A = tau_mag_A + tau_srp_A + tau_aero_A;
tau_total_B = tau_mag_B + tau_srp_B + tau_aero_B;

%% Angular velocity plots
%for validation of the t=46min spikes 
figure('Name', 'Angular Velocities', 'Position', [50,50,1000,700]);
subplot(2,1,1);
plot(tB/60, w_hist_A, 'LineWidth', 1.5); 
grid on;
xlabel('Time (min)'); 
ylabel('\omega (rad/s)');
title('Satellite A — Pointing angular velocities');
legend('\omega_x', '\omega_y', '\omega_z');
subplot(2,1,2);
plot(tAtt_B/60, w_hist_B, 'LineWidth', 1.5); 
grid on;
xlabel('Time (min)'); 
ylabel('\omega (rad/s)');
title('Satellite B — Tumbling angular velocities (disturbance-driven)');
legend('\omega_x', '\omega_y', '\omega_z');

%% Sat A torque plots
%for validation of the behaviour 
figure('Name', 'Satellite A Torques', 'Position', [100,100,1400,900]);
lbl = {'\tau_x', '\tau_y', '\tau_z'};
subplot(3,2,1);
plot(tB/60, tau_mag_A*1e6,  'LineWidth',1.5); 
grid on;
xlabel('Time (min)'); 
ylabel('\tau (\muN\cdotm)'); 
title('Sat A — Magnetic'); 
legend(lbl);
subplot(3,2,3); 
plot(tB/60, tau_srp_A*1e6,  'LineWidth',1.5); grid on;

xlabel('Time (min)'); 
ylabel('\tau (\muN\cdotm)'); 
title('Sat A — SRP'); legend(lbl);
subplot(3,2,5);
plot(tB/60, tau_aero_A*1e6, 'LineWidth',1.5); 
grid on;
xlabel('Time (min)'); 
ylabel('\tau (\muN\cdotm)'); 
title('Sat A — Aero drag'); 
legend(lbl);
mag_A_mag  = vecnorm(tau_mag_A,  2,2);
mag_A_srp  = vecnorm(tau_srp_A,  2,2);
mag_A_aero = vecnorm(tau_aero_A, 2,2);
mag_A_tot  = vecnorm(tau_total_A,2,2);
subplot(3,2,[2 4 6]);
semilogy(tB/60, mag_A_mag*1e6,  'b', ...
         tB/60, mag_A_srp*1e6,  'g', ...
         tB/60, mag_A_aero*1e6, 'r', ...
         tB/60, mag_A_tot*1e6,  'k--', 'LineWidth', 1.5);
grid on; 
xlabel('Time (min)'); 
ylabel('|\tau| (\muN\cdotm) — log');
title('Sat A — Torque magnitudes');
legend('Magnetic', 'SRP', 'Aero drag', 'Total');

%% Sat B torque plots
%for validation of the behaviour 
figure('Name', 'Satellite B Torques', 'Position', [150,150,1400,900]);
subplot(3,2,1); plot(tAtt_B/60, tau_mag_B*1e6,  'LineWidth',1.5); grid on;
xlabel('Time (min)'); ylabel('\tau (\muN\cdotm)'); title('Sat B — Magnetic'); legend(lbl);
subplot(3,2,3); plot(tAtt_B/60, tau_srp_B*1e6,  'LineWidth',1.5); grid on;
xlabel('Time (min)'); ylabel('\tau (\muN\cdotm)'); title('Sat B — SRP'); legend(lbl);
subplot(3,2,5); plot(tAtt_B/60, tau_aero_B*1e6, 'LineWidth',1.5); grid on;
xlabel('Time (min)'); ylabel('\tau (\muN\cdotm)'); title('Sat B — Aero drag'); legend(lbl);
mag_B_mag = vecnorm(tau_mag_B,  2,2);
mag_B_srp = vecnorm(tau_srp_B,  2,2);
mag_B_aero = vecnorm(tau_aero_B, 2,2);
mag_B_tot = vecnorm(tau_total_B,2,2);
subplot(3,2,[2 4 6]);
semilogy(tAtt_B/60, mag_B_mag*1e6,  'b', ...
         tAtt_B/60, mag_B_srp*1e6,  'g', ...
         tAtt_B/60, mag_B_aero*1e6, 'r', ...
         tAtt_B/60, mag_B_tot*1e6,  'k--', 'LineWidth', 1.5);
grid on; xlabel('Time (min)'); ylabel('|\tau| (\muN\cdotm) — log');
title('Sat B — Torque magnitudes');
legend('Magnetic', 'SRP', 'Aero drag', 'Total');


%% Print summary for reference

fprintf('TORQUE SUMMARY  (muN.m = 10^-6 N.m)\n');
fprintf('\nSatellite A (pointing chaser, 56 kg box):\n');
fprintf('  Magnetic: avg %7.4f  peak %7.4f  muN.m\n', mean(mag_A_mag)*1e6, max(mag_A_mag)*1e6);
fprintf('  SRP: avg %7.4f  peak %7.4f  muN.m\n', mean(mag_A_srp)*1e6, max(mag_A_srp)*1e6);
fprintf('  Aero drag: avg %7.4f  peak %7.4f  muN.m\n', mean(mag_A_aero)*1e6, max(mag_A_aero)*1e6);
fprintf('  Total: avg %7.4f  peak %7.4f  muN.m\n', mean(mag_A_tot)*1e6,  max(mag_A_tot)*1e6);
fprintf('\nSatellite B (tumbling debris, 45 kg rocket body):\n');
fprintf('  Magnetic: avg %7.4f  peak %7.4f  muN.m\n', mean(mag_B_mag)*1e6, max(mag_B_mag)*1e6);
fprintf('  SRP: avg %7.4f  peak %7.4f  muN.m\n', mean(mag_B_srp)*1e6, max(mag_B_srp)*1e6);
fprintf('  Aero drag: avg %7.4f  peak %7.4f  muN.m\n', mean(mag_B_aero)*1e6, max(mag_B_aero)*1e6);
fprintf('  Total: avg %7.4f  peak %7.4f  muN.m\n', mean(mag_B_tot)*1e6,  max(mag_B_tot)*1e6);

%% Distance between satellites over time
range_AB = vecnorm(sA_on_tB(:,1:3) - [xB, yB, zB], 2, 2);
figure('Name', 'Inter-satellite Distance', 'Position', [100, 100, 1000, 450]);
plot(tB/60, range_AB/1e3, 'c', 'LineWidth', 2);
grid on;
xlabel('Time (min)'); ylabel('Distance (km)');
title('Distance between Global-12 and Electron Kick Stage over time');
[min_range, idx_min] = min(range_AB);
[max_range, idx_max] = max(range_AB);
hold on;
plot(tB(idx_min)/60, min_range/1e3, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
plot(tB(idx_max)/60, max_range/1e3, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
legend('Inter-satellite range', ...
       sprintf('Min: %.1f km at t=%.1f min', min_range/1e3, tB(idx_min)/60), ...
       sprintf('Max: %.1f km at t=%.1f min', max_range/1e3, tB(idx_max)/60), ...
       'Location', 'best');
fprintf('\nInter-satellite distance:\n');
fprintf('  Min: %.2f km at t = %.2f min\n', min_range/1e3, tB(idx_min)/60);
fprintf('  Max: %.2f km at t = %.2f min\n', max_range/1e3, tB(idx_max)/60);
fprintf('  Mean: %.2f km\n', mean(range_AB)/1e3);


%% Orbital Energy check
%(Change number of orbits to >1 to see effects)
r_mag = vecnorm(sA(:,1:3), 2, 2);
v_mag = vecnorm(sA(:,4:6), 2, 2);
energy = 0.5*v_mag.^2 - mu./r_mag;
figure; plot(tA/60, energy); grid on;
xlabel('Time (min)'); ylabel('Specific energy (J/kg)');
title('Orbital energy'); % should oscilate due to the J2 perterbations adding and removing enrgy during the orbit 
%and gradually decrease over time due to drag removing energy,


%%  FUNCTIONS
% Orbital dynamics — two-body gravity + J2 + atmospheric drag + SRP
function dstatedt = OrbitsGeneral(t, s, sat, J2, Re, mu)
    r = s(1:3);   v = s(4:6);
    x = r(1); y = r(2); z = r(3);   rn = norm(r);

    % Point mass gravity
    a_2body = -mu * r / rn^3;

    % J2 — Earth oblateness, (
    % Causes RAAN drift (~7 deg/day at 400 km) and argument of perigee precession
    fJ2 = (3/2) * J2 * mu * Re^2 / rn^7;
    a_J2 = fJ2 * [x*(x^2+y^2-4*z^2); y*(x^2+y^2-4*z^2); z*(3*x^2+3*y^2-2*z^2)];

    % Atmospheric drag — inertial velocity used directly (no co-rotation)
    v_norm   = norm(v);
    altitude = rn - 6.371e6;
    rho = atmosphericDensity(altitude);
    a_drag = -0.5 * sat.Cd * sat.A_drag * rho * v_norm * v / sat.m;

    % Solar radiation pressure
    r_sun = getSunPosition(t);
    a_srp = solarRadiationPressure(r, r_sun, sat.Cr, sat.A_srp, sat.m);

    dstatedt = [v; a_2body + a_J2 + a_drag + a_srp];
end


% Piecewise exponential atmosphere — based on CIRA-2012 table values
function rho = atmosphericDensity(alt_m)
    layers = [
        200e3, 2.840e-10, 35237.95;
        220e3, 1.610e-10, 38680.52;
        240e3, 9.600e-11, 42103.83;
        260e3, 5.970e-11, 45057.01;
        280e3, 3.830e-11, 43563.75;
        300e3, 2.420e-11, 55704.25;
        320e3, 1.690e-11, 53147.88;
        340e3, 1.160e-11, 53646.01;
        360e3, 7.990e-12, 56270.79;
        380e3, 5.600e-12, 57716.30;
        400e3, 3.960e-12, 59529.60;
        420e3, 2.830e-12, 60197.28;
        440e3, 2.030e-12, 61962.98;
        460e3, 1.470e-12, 62971.55;
        480e3, 1.070e-12, 64572.33;
        500e3, 7.850e-13, 65336.02;
        520e3, 5.780e-13, 67087.77;
        540e3, 4.290e-13, 67506.94;
        560e3, 3.190e-13, 69269.45;
        580e3, 2.390e-13, 70545.07;
        600e3, 1.800e-13, 71351.62;
        620e3, 1.360e-13, 74553.43;
        640e3, 1.040e-13, 75509.48;
        660e3, 7.980e-14, 77261.35;
        680e3, 6.160e-14, 80172.90;
        700e3, 4.800e-14, 81901.10;
        720e3, 3.760e-14, 86023.11;
        740e3, 2.980e-14, 88958.94;
        760e3, 2.380e-14, 93120.58;
        780e3, 1.920e-14, 99379.10;
        800e3, 1.570e-14, 99379.10];
    idx = find(alt_m >= layers(:,1), 1, 'last');
% if altitude is below the lowest layer just use the bottom value
if isempty(idx)
    idx = 1; 
end
    if idx > size(layers,1)
        idx = size(layers,1); 
    end
    rho = layers(idx,2) * exp(-(alt_m - layers(idx,1)) / layers(idx,3));
    rho = max(rho, 1e-15);
end

% Simplified sun position — circular orbit with 23.4 deg ecliptic obliquity
function r_sun = getSunPosition(t)
    AU = 1.496e11;
    w_earth = 2*pi / (365.25*24*3600);
    theta = w_earth * t;
    eps = deg2rad(23.4);
    r_sun = AU * [cos(theta); sin(theta)*cos(eps); sin(theta)*sin(eps)];
end

% SRP translational acceleration with dual-cone shadow check
function a_srp = solarRadiationPressure(r_sat, r_sun, Cr, A_srp, m_sat)
    AU  = 1.496e11;   P_sun = 4.56e-6;
    rs = r_sun - r_sat;   ds = norm(rs);
    nu = shadowFunction(r_sat, r_sun, 6.371e6);
    a_srp = nu * Cr * (A_srp/m_sat) * P_sun * (AU/ds)^2 * (rs/ds);
end

% Dual-cone shadow model
% Returns nu = 1 (sunlit), 0 (umbra), or 0 < nu < 1 (penumbra)
%should expect roughly 34min of cutoff with no SRP on graphs
function nu = shadowFunction(r_sat, r_sun, R_earth)
    R_sun = 6.957e8;
    d_sat = norm(r_sat);
    d_ss = norm(r_sun - r_sat);
    a_sun  = asin(min(1, R_sun   / d_ss));
    a_Earth = asin(min(1, R_earth / d_sat));
    cos_sep = dot(-r_sat, r_sun - r_sat) / (d_sat * d_ss);
    sep = acos(max(-1, min(1, cos_sep)));
    if sep >= a_sun + a_Earth
        nu = 1;
    elseif sep <= abs(a_Earth - a_sun)
        nu = 0;
    else
        nu = (sep - abs(a_Earth-a_sun)) / (a_sun+a_Earth - abs(a_Earth-a_sun));
    end
end

% Earth magnetic field — tilted dipole (IGRF-13, 11.5 deg offset)
function B = earthMagneticField(r_sat)
    mu_0  = 4*pi*1e-7;
    M_dip = 7.94e22;
    tilt = deg2rad(11.5);
    m_hat = [sin(tilt); 0; cos(tilt)];
    rn = norm(r_sat);
    r_hat = r_sat / rn;
    coeff = (mu_0/(4*pi)) * M_dip / rn^3;
    B = coeff * (3*dot(m_hat, r_hat)*r_hat - m_hat);
end

% SRP torque 
function tau_srp = srpTorque(r_sat, q, r_sun, sat)
    AU = 1.496e11;   P_sun = 4.56e-6;
    rs = r_sun - r_sat;   ds = norm(rs);   s_I = rs/ds;
    nu  = shadowFunction(r_sat, r_sun, 6.371e6);
    if nu == 0, tau_srp = [0;0;0]; 
        return; 
    end
    C = quat2dcm(q);
    s_B = C' * s_I;
    P = nu * P_sun * (AU/ds)^2;
    tau_srp = [0;0;0];
    for i = 1:length(sat.faces)
        n_i = sat.faces{i};
        ct  = dot(n_i, s_B);
        if ct > 0
            rho_i = sat.reflectivity(i) - 1;
            F_i = P * sat.areas(i) * ct * ((1-rho_i)*s_B + 2*rho_i*ct*n_i);
            tau_srp = tau_srp + cross(sat.cp_srp{i}, F_i);
        end
    end
end

% Aerodynamic drag torque — dynamic pressure on each opposing face
function tau_aero = aeroDragTorque(r_sat, q, v_inertial, sat)
    altitude = norm(r_sat) - 6.371e6;
    rho = atmosphericDensity(altitude);
    vn = norm(v_inertial);
    if vn < 1e-3, tau_aero = [0;0;0]; 
        return;
    end
    C = quat2dcm(q);
    v_hat_B = C' * (v_inertial / vn);
    q_dyn = 0.5 * sat.Cd * rho * vn^2;
    tau_aero = [0;0;0];
    for i = 1:length(sat.faces)
        n_i = sat.faces{i};
        ct  = -dot(n_i, v_hat_B);
        if ct > 0
            F_i = -q_dyn * sat.areas(i) * ct * v_hat_B;
            tau_aero = tau_aero + cross(sat.cp_aero{i}, F_i);
        end
    end
end

% Attitude dynamics — magnetic + SRP + aero torques drive Sat B tumbling
function dstate = AttitudeDynamics(t, state, orb_self, sat)
    q = state(1:4);  
    w = state(5:7);
    q_err = norm(q)^2 - 1;
    q = q / norm(q);
    sv = orb_self(t);
    r_sat = sv(1:3);  
    v_sat = sv(4:6);
    I = sat.I;
    B_I = earthMagneticField(r_sat);
    C = quat2dcm(q);
    tau_mag = cross(sat.mag_diop, C'*B_I);
    tau_srp = srpTorque(r_sat, q, getSunPosition(t), sat);
    tau_aero = aeroDragTorque(r_sat, q, v_sat, sat);
    tau_total = tau_mag + tau_srp + tau_aero;
    % Quaternion kinematics with Baumgarte stabilisation (k = 1 s^-1)
    Omega = [ 0 , -w(1), -w(2), -w(3);
            w(1), 0, w(3), -w(2);
            w(2), -w(3), 0, w(1);
            w(3), w(2), -w(1), 0 ];
    qdot  = 0.5 * Omega * q - 1.0 * q_err * q;
    % Euler's rotational equations
    wdot  = I \ (tau_total - cross(w, I*w));
    dstate = [qdot; wdot];
end

% Quaternion to DCM — scalar-first q = [q0, q1, q2, q3]
function C = quat2dcm(q)
    q  = q / norm(q);
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);
    C = [1-2*(q2^2+q3^2), 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2);
         2*(q1*q2+q0*q3), 1-2*(q1^2+q3^2), 2*(q2*q3-q0*q1);
         2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 1-2*(q1^2+q2^2)];
end

% DCM to quaternion — Shepperd's method, NaN/Inf guarded
function q = dcm2quat(C)

    tr = trace(C);
    if tr > 0
        s  = 0.5 / sqrt(tr + 1);   
        q0 = 0.25 / s;
        q1 = (C(3,2)-C(2,3))*s;   
        q2 = (C(1,3)-C(3,1))*s;   
        q3 = (C(2,1)-C(1,2))*s;
    elseif C(1,1) > C(2,2) && C(1,1) > C(3,3)
        s  = 2 * sqrt(max(1e-12, 1+C(1,1)-C(2,2)-C(3,3)));
        q0 = (C(3,2)-C(2,3))/s; 
        q1 = 0.25*s;
        q2 = (C(1,2)+C(2,1))/s; 
        q3 = (C(1,3)+C(3,1))/s;
    elseif C(2,2) > C(3,3)
        s  = 2 * sqrt(max(1e-12, 1+C(2,2)-C(1,1)-C(3,3)));
        q0 = (C(1,3)-C(3,1))/s;  
        q1 = (C(1,2)+C(2,1))/s;
        q2 = 0.25*s; 
        q3 = (C(2,3)+C(3,2))/s;
    else
        s  = 2 * sqrt(max(1e-12, 1+C(3,3)-C(1,1)-C(2,2)));
        q0 = (C(2,1)-C(1,2))/s; 
        q1 = (C(1,3)+C(3,1))/s;
        q2 = (C(2,3)+C(3,2))/s; 
        q3 = 0.25*s;
    end
    q = [q0; q1; q2; q3];
    q = q / norm(q);
end

%% convert to JSON file
simData = struct();
simData.satelliteA.orbit.time = tA;
simData.satelliteA.orbit.position.xA = xA;
simData.satelliteA.orbit.position.yA  = yA;
simData.satelliteA.orbit.position.zA = zA;
simData.satelliteA.orbit.velocity.vxA = sA(:,4);
simData.satelliteA.orbit.velocity.vyA = sA(:,5);
simData.satelliteA.orbit.velocity.vzA = sA(:,6);
simData.satelliteA.attitude.quaternion.q0 = q_hist_A(:,1);
simData.satelliteA.attitude.quaternion.q1 = q_hist_A(:,2);
simData.satelliteA.attitude.quaternion.q2 = q_hist_A(:,3);
simData.satelliteA.attitude.quaternion.q3 = q_hist_A(:,4);
 
simData.satelliteB.orbit.time  = tB;
simData.satelliteB.orbit.position.x = xB;
simData.satelliteB.orbit.position.y = yB;
simData.satelliteB.orbit.position.z = zB;
simData.satelliteB.orbit.velocity.vx = sB(:,4);
simData.satelliteB.orbit.velocity.vy = sB(:,5);
simData.satelliteB.orbit.velocity.vz = sB(:,6);
simData.satelliteB.attitude.quaternion.q0 = q_hist_B(:,1);
simData.satelliteB.attitude.quaternion.q1 = q_hist_B(:,2);
simData.satelliteB.attitude.quaternion.q2  = q_hist_B(:,3);
simData.satelliteB.attitude.quaternion.q3 = q_hist_B(:,4);
simData.satelliteB.attitude.angularVelocity.wx = w_hist_B(:,1);
simData.satelliteB.attitude.angularVelocity.wy = w_hist_B(:,2);
simData.satelliteB.attitude.angularVelocity.wz = w_hist_B(:,3);
 
simData.parameters.EarthRadius = R;
simData.parameters.earthMass = M_earth;         
simData.parameters.gravitationalConstant = G;
simData.parameters.mu  = mu;
simData.parameters.orbitA.semiMajorAxis = semi_majorA;
simData.parameters.orbitA.eccentricity = eccenA;
simData.parameters.orbitA.inclination  = incA;
simData.parameters.orbitB.semiMajorAxis = semi_majorB;
simData.parameters.orbitB.eccentricity = eccenB;
simData.parameters.orbitB.inclination  = incB;
simData.parameters.InitialAltitude= Altitude_ref;    
 
simData.metadata.simulationTime = datestr(now);
simData.metadata.numberOfTimeSteps = length(tB);
simData.metadata.orbitPeriodA = Ta;
simData.metadata.orbitPeriodB = Tb;
 
jsonStr = jsonencode(simData);
filename = 'orbital_simulation_data.json';
fid = fopen(filename, 'w');
fprintf(fid, '%s', jsonStr);
fclose(fid);
fprintf('Simulation data exported to: %s\n', filename);

%TO DO: Atmosphere model only goes up to 800km, works for LEO and
%simulations but might break using values for apogee and perigee above 800
% Add J3 and other perturbation values for more accurate model
% Fix drag coefficient value of 2.2 (Used as and approximation but likely
% off by 5-10% of real value