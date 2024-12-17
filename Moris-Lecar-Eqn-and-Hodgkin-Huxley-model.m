clear; clc;

global C;
global gCa;
global VCa;
global gK;
global VK;
global gL;
global VL;
global v1;
global v2;
global v3;
global v4;
global phi;
global Iext;

gCa = 4.4;
gK = 8.0;
gL = 2.0;
VCa = 120.0;
VK = -84.0;
VL = -60.0;
v1 = -1.2;
v2 = 18.0;
v3 = 2.0;
v4 = 30.0;
phi = 0.02;
C = 20.0;
Iext = 0;

% --------------------- QUESTION 2 --------------------------------

% [2] Nullclines for Moris Lecar equation

% m_inf = 0.5.*(1 + tanh((V-v1)./v2))
V = linspace(-75, 80);
V_nullcline = (Iext - gCa.*(0.5.*(1 + tanh((V-v1)./v2))).*(V-VCa) - gL.*(V-VL))./(gK.*(V-VK));
w_nullcline = (0.5.*(1 + tanh((V-v3)./v4)));

figure;
plot(V, V_nullcline);
hold on;
plot(V, w_nullcline);
hold on;

xlabel('V'); 
ylabel('w');
title('(2) Phase Plane Plot of V and w nullclines');
grid on;

% [2] Finding equilibrium points for the MLE

m_inf = @(V) 0.5 * (1 + tanh((V - v1) / v2)); % m-infinity
w_nullcline = @(V) 0.5 * (1 + tanh((V - v3) / v4)); % w-nullcline

F = @(V) (1 / C) * (Iext - gCa * m_inf(V) * (V - VCa) ...
        - gK * w_nullcline(V) * (V - VK) - gL * (V - VL));

V_eq = fzero(F, [-80, 60]);

w_eq = w_nullcline(V_eq);

fprintf("----- Question 1 -------- \n");
disp('Equilibrium point (using fzero):');
disp(['V_eq = ', num2str(V_eq)]);
disp(['w_eq = ', num2str(w_eq)]);


syms V w
Vnc_eqn = (1/C)*(Iext - gCa*(0.5*(1 + tanh((V - v1)/v2)))*(V - VCa) - gK*w*(V - VK) - gL*(V - VL)) == 0;
wnc_eqn = (0.5*(1 + tanh((V - v3)/v4)) - w) == 0;
eq_pt_0 = vpasolve([Vnc_eqn, wnc_eqn], [V, w]);

V_eq = double(eq_pt_0.V);
w_eq = double(eq_pt_0.w);

fprintf('Equilibrium point (using vpasolve): \n');
fprintf('V = %.4f\n', V_eq);
fprintf('w = %.4f\n', w_eq);


% [2] ------ Quiver plot -------
x = linspace(-80, 80, 16);
y = linspace(0, 1, 12);

[V_quiv, w_quiv] = meshgrid(x, y);

dV_dt = (1/C) * (Iext + gCa * (0.5 * (1 + tanh((V_quiv - v1) / v2))) .* (VCa - V_quiv) + gK * w_quiv .* (VK - V_quiv) + gL * (VL - V_quiv));
dw_dt = phi * ((0.5 * (1 + tanh((V_quiv - v3) / v4))) - w_quiv) .* cosh((V_quiv - v3) / (2 * v4));

quiver(V_quiv, w_quiv, dV_dt, dw_dt);
axis tight;


% --------------------- QUESTION 3 ----------------------------------------
% [3] --------- Jacobian and Eigenvalues of the equilibrium point ---------

syms V w
f1 = (1/C).*(Iext - gCa .*(0.5 .*(1 + tanh((V - v1) ./ v2))) .* (V - VCa) - gK .* w .* (V - VK) - gL .* (V - VL));
f2 = phi * ((0.5 .* (1 + tanh((V - v3) ./ v4))) - w) .* cosh((V - v3) ./ (2 .* v4));

% JACOBIAN
J = jacobian([f1; f2], [V, w]);

% Assuming V_eq and w_eq are already defined as equilibrium points
J_equi = subs(J, {V, w}, {V_eq, w_eq});  % Calculate Jacobian at equilibrium points
J_equi = double(J_equi);

% Calculate eigenvalues
eigenvalues = eig(J_equi);

lambda1 = eigenvalues(1);
lambda2 = eigenvalues(2);

% Print Jacobian and Eigenvalues
fprintf("\n ----- Question 2 -------- \n")

fprintf("Jacobian @ [%.4f, %.4f]: \n", V_eq, w_eq);  % Print equilibrium points
disp(J_equi);
disp("Eigenvalues: ")
disp(eigenvalues);

equilibriumType = type_of_equilibrium(lambda1, lambda2);

% Print the result
fprintf('The type of equilibrium is: %s\n', equilibriumType);

% ----------------------- DONE -----------------------------
% ----------------------------------------------------------

% ---------------- QUESTION 5 -----------------------------

function dSdt = morris_lecar(t, S)
    global C gCa VCa gK VK gL VL v1 v2 v3 v4 phi Iext;

    V = S(1); % Membrane potential
    w = S(2); % Gating variable

    % Define the Morris-Lecar equations
    dV_dt = (1/C) * (Iext - gCa * (0.5 * (1 + tanh((V - v1) / v2))) * (V - VCa) ...
                     - gK * w * (V - VK) - gL * (V - VL));
    dw_dt = phi * (0.5 * (1 + tanh((V - v3) / v4)) - w) * cosh((V - v3) / (2 * v4));

    % Return the derivatives as a column vector
    dSdt = [dV_dt; dw_dt];

end

% Define the Morris-Lecar equations as a function
options = odeset('RelTol',1e-3,'AbsTol',1e-6, 'refine',5, 'MaxStep', 1);

Iext = 0;
tSpan = [0, 300];
initial = [0,w_eq];

phi = 0.02;
[t1, S1] = ode15s(@(t,S)morris_lecar(t,S),tSpan, initial, options);

phi = 0.04;
[t2, S2] = ode15s(@(t,S)morris_lecar(t,S),tSpan, initial, options);

phi = 0.01;
[t3, S3] = ode15s(@(t,S)morris_lecar(t,S),tSpan, initial, options);

phi = 0.02;

% Plot the Action Potentials --------------------------------------
figure;
hold on;
plot(t1,S1(:,1));
plot(t2,S2(:,1));
plot(t3,S3(:,1));
xlabel('Time(in ms)');
ylabel('Volatage(in mV)');
title('(5a) Action potentials with different \phi');
legend('\phi = 0.01','\phi = 0.02','\phi = 0.04');
grid on;

% Plot the PHASE PLANE PLOT -----------------------------------------
figure;
hold on
Vnc = @(V) (Iext - gCa .* (0.5 .* (1 + tanh((V-v1)./v2))) .* (V-VCa) - gL .* (V-VL))./(gK .* (V-VK));
wnc = @(V) (0.5 .* (1 + tanh((V-v3)./v4)));
fplot(@(V) Vnc(V), [-80 100],'k');
fplot(@(V) wnc(V), [-80 100],'k');

plot(S1(:,1),S1(:,2));
plot(S2(:,1),S2(:,2));
plot(S3(:,1),S3(:,2));
xlabel('V(in mV)');
ylabel('w');
ylim([0,1]);
title('(5b) Phase Plane Plot (MLE)');
legend('V nullcline','w nullcline','\phi = 0.01','\phi = 0.02','\phi = 0.04');
grid on;

% --------------- DONE ------------------------------------------


%--------------- QUESTION 6 -------------------------------
% Simulate depolarizing current pulses of various amplitudes by setting
% Voltages to relatively positive values with starting "w" at equilibrium

Iext = 0.0;
tSpan = [0 400]; 
V_initial = linspace(-60, 0, 400);
max_amplitude = zeros(1, 400);
point = 0;

for i = 1:400
    [t, S] = ode15s(@(t,S)morris_lecar(t,S),tSpan, [V_initial(i),w_eq], options);
    max_V_value = max(S(:,1));
    max_amplitude(i) = max_V_value;
    if(max_amplitude(i) >= 0)
        threshold = V_initial(i);
        point = i;
        break;
    end
end

fprintf("\n ------ Question 6 ------------------- \n");
fprintf("Threshold potential under Depolarising current = %.2f mV \n", threshold);

figure;
plot(V_initial, max_amplitude, 'r', 'LineWidth', 0.5);
hold on;
xline(threshold, '-b', 'LineWidth', 1, 'Label', 'Threshold', 'LabelVerticalAlignment', 'bottom', 'LabelHorizontalAlignment', 'center');
xlabel('Initial Voltage');
ylabel('Amplitude of voltage during depolarization');
title('(6) Depolarizing current pulses');
grid on;

% ------------------- QUESTION 7 --------------------------------
% ---------------------------------------------------------------

% CREATE FUNCTION

function equilibriumType = type_of_equilibrium(lambda1, lambda2)
    % Check the TYPE OF EQUILIBRIUM
    if isreal(lambda1) && isreal(lambda2) % Both eigenvalues are real
        if (lambda1 < 0) && (lambda2 < 0)
            equilibriumType = 'Stable Equilibrium';
        elseif (lambda1 > 0) && (lambda2 > 0)
            equilibriumType = 'Unstable Equilibrium';
        else
            equilibriumType = 'Saddle Point';
        end
    
    else % Complex eigenvalues
        if real(lambda1) < 0 && real(lambda2) < 0
            equilibriumType = 'Stable Spiral';
        elseif real(lambda1) > 0 && real(lambda2) > 0
            equilibriumType = 'Unstable Spiral';
        else
            equilibriumType = 'Saddle Focus';
        end
    end
end


% Run the model with I_ext = 86 μA/cm2 with three sets of initial conditions:
% (i) Initial conditions = equilibrium conditions @ I_ext = 0.0
% (ii) Initial conditions = equilibrium conditions @ I_ext = 86.0 μA/cm2
% (iii) Initial conditions = (-27.9, 0.17)

tSpan = [0 600];

fprintf("\n ---------- QUESTION 7 -------------------------- \n");

% (i) Initial conditions = equilibrium conditions @ I_ext = 0.0
fprintf("For Iext = 0 \n");
fprintf("Equilibrium points: \n");
fprintf("V = %.4f mV \n", V_eq);
fprintf("w = %.4f \n", w_eq);


% -----------------------------------------------------------------


% (ii) Initial conditions = equilibrium conditions @ I_ext = 86.0 μA/cm2
Iext = 86.0;
syms V w
Vnc_eqn = (1/C) .* (Iext - gCa .* (0.5 .* (1 + tanh((V - v1)/v2))) .* (V - VCa) - gK .* w .*(V - VK) - gL .* (V - VL)) == 0;
wnc_eqn = (0.5 .* (1 + tanh((V - v3)./v4)) - w) == 0;
eq_pt_0 = vpasolve([Vnc_eqn, wnc_eqn], [V, w]);

V_eq_2 = double(eq_pt_0.V);
w_eq_2 = double(eq_pt_0.w);

fprintf("\nFor Iext = 86.0 μA/cm2 \n");
fprintf("Equilibrium points: \n");
fprintf("V = %.4f mV \n", V_eq_2);
fprintf("w = %.4f \n", w_eq_2);

initial = [V_eq_2, w_eq_2];
[t2, S2] = ode15s(@(t,S)morris_lecar(t,S), tSpan, initial, options);

% Characterize equilibrium points for Iext = 86.0 -------------------------
Iext = 86.0;
syms V w
f1 = (1/C)*(Iext - gCa*(0.5*(1 + tanh((V - v1) / v2))) * (V - VCa) - gK * w * (V - VK) - gL * (V - VL));
f2 = phi * ((0.5 * (1 + tanh((V - v3) / v4))) - w) * cosh((V - v3) / (2 * v4));

% JACOBIAN
J = jacobian([f1; f2], [V, w]);

% Assuming V_eq and w_eq are already defined as equilibrium points
J_equi = subs(J, {V, w}, {V_eq_2, w_eq_2});  % Calculate Jacobian at equilibrium points
J_equi = double(J_equi);

% Calculate eigenvalues
eigenvalues = eig(J_equi);

lambda1 = eigenvalues(1);
lambda2 = eigenvalues(2);

fprintf("Jacobian @ [%.4f, %.4f]: \n", V_eq_2, w_eq_2);  % Print equilibrium points
disp(J_equi);
disp("Eigenvalues: ")
disp(eigenvalues);

equilibriumType = type_of_equilibrium(lambda1, lambda2);

% Print the result
fprintf('The type of equilibrium is: %s \n', equilibriumType);

% PHASE PLANE PLOTS ------------------------------------------------
Iext = 86.0;
tSpan = [0, 600];

initial = [V_eq, w_eq];
[t1, S1] = ode15s(@(t,S)morris_lecar(t,S), tSpan, initial, options);

initial = [V_eq_2 + 0.1, w_eq_2 + 0.001];
[t2, S2] = ode15s(@(t,S)morris_lecar(t,S), tSpan, initial, options);

initial = [-27.9, 0.17];
[t3, S3] = ode15s(@(t,S)morris_lecar(t,S), tSpan, initial, options);

figure;
hold on;

% Nullclines for Iext = 0
Iext_0 = 0.0; 
Vnc_0 = @(V) (Iext_0 - gCa.*(0.5.*(1+tanh((V-v1)./v2))).*(V-VCa) - gL.*(V-VL))./(gK.*(V-VK));
fplot(@(V) Vnc_0(V), [-80 100],'k--');  % Dashed line for distinction

% Nullclines for Iext = 86
Iext_86 = 86.0;
Vnc_86 = @(V) (Iext_86 - gCa.*(0.5.*(1+tanh((V-v1)./v2))).*(V-VCa) - gL.*(V-VL))./(gK.*(V-VK));
fplot(@(V) Vnc_86(V), [-80 100],'k');

wnc = @(V) (0.5.*(1+tanh((V-v3)./v4)));
fplot(@(V) wnc(V), [-80 100],'m');

% Plot trajectories for each condition
plot(S1(:,1), S1(:,2), 'r', 'LineWidth', 1.5); % Iext = 0 equilibrium point
plot(S2(:,1), S2(:,2), 'g', 'LineWidth', 1.5); % Iext = 86 equilibrium point
plot(S3(:,1), S3(:,2), 'b', 'LineWidth', 1.5); % Iext = 86 random initial point
plot(-60.8554, 0.0149, 'bo', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
plot(-27.9, 0.17, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6);

xlim([-80 100]);
ylim([0 1]);
xlabel('V (mV)');
ylabel('w');
title('(7) Phase plane plot for different conditions');
legend('V nullcline I_{ext} = 0','V nullcline I_{ext} = 86','w nullcline', ...
       'Equilibrium point for I_{ext} = 0','Equilibrium point for I_{ext} = 86 ', ...
       'Random point');
grid on;

% ------------------------- QUESTION 8 ------------------------------
% ------------------------------------------------------------------
% Find the contour that divides the phase plane into those initial conditions that 
% converge to the equilibrium point and those that converge to the limit cycle

fprintf("\n --------------- QUESTION 8 -------------------- \n");

Iext = 86;

% Plotting the Phase Plane with V and w null clines
figure;
hold on;

Vnc1 = @(V) (Iext - gCa .* (0.5.*(1 + tanh((V-v1)./v2))) .* (V-VCa) - gL .* (V-VL))./(gK .* (V-VK));
wnc1 = @(V) (0.5 .* (1 + tanh((V-v3)./v4)));
fplot(@(V) Vnc1(V), [-80 100], 'k');
fplot(@(V) wnc1(V), [-80 100], 'k');

% Equilibrium point for the system
syms V w
Vnc1_eqn = (1/C)*(Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gK*w*(V-VK) - gL*(V-VL)) == 0;
wnc1_eqn = (0.5*(1+tanh((V-v3)/v4)) - w) == 0;
eq_pt_1 = vpasolve([Vnc1_eqn, wnc1_eqn], [V, w]);

V_eq1 = double(eq_pt_1.V);
w_eq1 = double(eq_pt_1.w);

plot(V_eq1, w_eq1, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 6);
xlabel('V (mV)');
ylabel('w');
fprintf("Equilibrium point: \n");
fprintf("V = %.4f mV \n", V_eq1); 
fprintf("w = %.4f \n", w_eq1);

% Run the model backwards in time
tSpan1 = [0, 1000];
tSpan2 = [0, -600];

[t1,S1] = ode45(@(t,S)morris_lecar(t,S), tSpan1, [V_eq, w_eq]);
[t2,S2] = ode45(@(t,S)morris_lecar(t,S), tSpan1, [-30, 0.3]);
[t3,S3] = ode45(@(t,S)morris_lecar(t,S), tSpan1, [V_eq1 + 0.1, w_eq1 + 0.001]);
[t4,S4] = ode45(@(t,S)morris_lecar(t,S), tSpan2, [-30, 0.3]);

plot(S1(:,1), S1(:,2), 'r', 'LineWidth', 1);
plot(S2(:,1), S2(:,2), 'g', 'LineWidth', 1);
plot(S3(:,1), S3(:,2), 'b', 'LineWidth', 1);
plot(S4(:,1), S4(:,2), 'm', 'LineWidth', 1);
xlim([-80 100]);
ylim([0 1]);
xlabel('V (mV)');
ylabel('w');
legend('V nullcline', 'w nullcline', 'Equilibrium point', 'Equi. points of I_{ext} = 0', ...
    'Random point (Limit cycle)', 'Equi. points of I_{ext} = 86 (Converge to equilibrium point)', 'Random point (backwards in time)');
title('(8) Simulation ran backwards in time');
grid on;

[t1,S1] = ode45(@(t,S)morris_lecar(t,S), tSpan1, [V_eq, w_eq]);
[t2,S2] = ode45(@(t,S)morris_lecar(t,S), tSpan1, [-30, 0.15]);
[t3,S3] = ode45(@(t,S)morris_lecar(t,S), tSpan1, [V_eq1 + 0.1, w_eq1 + 0.001]);
[t4,S4] = ode45(@(t,S)morris_lecar(t,S), tSpan2, [-30, 0.15]);

figure;
hold on;

fplot(@(V) Vnc1(V), [-80 100], 'k');
fplot(@(V) wnc1(V), [-80 100], 'k');
plot(S1(:,1), S1(:,2), 'r', 'LineWidth', 1);
plot(S2(:,1), S2(:,2), 'g', 'LineWidth', 1);
plot(S3(:,1), S3(:,2), 'b', 'LineWidth', 1);
plot(S4(:,1), S4(:,2), 'm', 'LineWidth', 1);

xlim([-80 100]);
ylim([0 1]);
xlabel('V (mV)');
ylabel('w');
legend('V nullcline', 'w nullcline', 'Equi. points of I_{ext} = 0', ...
    'Random point (Converge to equilibrium point)', 'Equi. points of I_{ext} = 86 (Converge to equilibrium point)', 'Random point (UPO backwards in time)');
title('(8) Simulation ran backwards in time 2');
grid on;

tSpan1 = [0, 1000];
[t1,S1] = ode45(@(t,S)morris_lecar(t,S), tSpan1, [V_eq, w_eq]);
[t2,S2] = ode45(@(t,S)morris_lecar(t,S), tSpan1, [-30, 0.15]);
figure;
plot(t1, S1(:,1), "Color",'r');
hold on;
plot(t2, S2(:,1), "Color",'g');
grid on;
xlabel('time (t)');
ylabel('Membrane Potential (mV)');
title('(8) Simulation of subthreshold vs action potential');
legend('Sustained Action Potential', 'Subthreshold Underdamped behaviour')

% -------------------- DONE --------------------------------------
% ----------------------------------------------------------------

% ------------------ QUESTION 9 ---------------------------------
% ---------------------------------------------------------------

% Rate of firing action potentials
function F = firing_rate(t, S)
    n = size(S);
    n = n(1);
    S = S(:, 1);
    
    spikes = 0;
    % When membrane potential crosses 0 --> Action potential has occured
    for i = 1:n
        if (S(i) > 0)
            spikes = 1;
        end
    end

    % If there are no spikes, terminate the function
    if (spikes == 0)
        F = 0;
        return
    end
    
    % Move ahead till we get a negative signal: 

    % The time point at which the membrane potential crosses ZERO potential
    % at the point of repolarisation (going backwards in time);

    while S(n) > 0 
        n = n - 1;
    end
    
    % Record first positive signal: second time-point
    while S(n) < 0 
        n = n -1 ;
    end

    t2 = t(n);

    % Recording the time point at which membrane potential crosses ZERO
    % again at the point of depolarisation

    while (S(n) > 0)
        n = n-1;
    end

    % Record second positive signal: first time point

    while S(n) < 0
        n = n-1;
    end

    t1 = t(n);
    
    % Time scale of action potential is in miliseconds (ms)
    % Frequency of firng = 1000/time (Hz)
    F = 1000 / (t2 - t1);
end

fprintf("\n ------------------------- QUESTION 9 ------------------------------ \n");
Iext_arr = [80, 86, 90];

figure;
hold on;

for i = 1:3
    Iext = Iext_arr(i);
    tSpan = [0 1000];

    % Finding the equilibrium point for this system using MATLAB
    syms V w
    Vnc1_eqn = (1/C)*(Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gK*w*(V-VK) - gL*(V-VL)) == 0;
    wnc1_eqn = (0.5*(1+tanh((V-v3)/v4)) - w) == 0;
    eq_pt_1 = vpasolve([Vnc1_eqn, wnc1_eqn], [V, w]);

    V_eq_4 = double(eq_pt_1.V);
    w_eq_4 = double(eq_pt_1.w);

    [t1, S1] = ode15s(@(t,S)morris_lecar(t,S), tSpan, [V_eq_4 + 0.5, w_eq_4 + 0.05]);

    subplot(2, 2, i);
    fplot(@(V) Vnc1(V), [-80 100], 'k');
    fplot(@(V) wnc1(V), [-80 100], 'k');
    plot(V_eq_4 + 0.5, w_eq_4 + 0.05, 'bo', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
    hold on;
    plot(V_eq_4, w_eq_4, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 3);
    hold on;
    plot(S1(:,1), S1(:,2), 'r', 'LineWidth', 1);
    xlabel('V (mV)');
    ylabel('w');
    title(sprintf('(9) Plot for Iext = %d μA/cm²', Iext));
    grid on;

    % Stability Analysis @ equilibrium point
    fprintf('Equilibrium point for Iext = %d μA/cm2 @ (%.2f,%.2f) \n', Iext, V_eq_4, w_eq_4);

    syms V w
    f1 = (1/C)*(Iext - gCa*(0.5*(1 + tanh((V - v1) / v2))) * (V - VCa) - gK * w * (V - VK) - gL * (V - VL));
    f2 = phi * ((0.5 * (1 + tanh((V - v3) / v4))) - w) * cosh((V - v3) / (2 * v4));
    
    % JACOBIAN
    J = jacobian([f1; f2], [V, w]);
   
    J_equi = subs(J, {V, w}, {V_eq_4, w_eq_4});  % Calculate Jacobian at equilibrium points
    J_equi = double(J_equi);
    
    % Calculate eigenvalues
    eigenvalues = eig(J_equi);
    
    lambda1 = eigenvalues(1);
    lambda2 = eigenvalues(2);
    
    fprintf("Jacobian @ [%.4f, %.4f]: \n", V_eq_4, w_eq_4);  % Print equilibrium points
    disp(J_equi);
    disp("Eigenvalues: ");
    disp(eigenvalues);
    
    % Check the TYPE OF EQUILIBRIUM
    if isreal(lambda1) && isreal(lambda2) % Both eigenvalues are real
        if (lambda1 < 0) && (lambda2 < 0)
            equilibriumType = 'Stable Equilibrium';
        elseif (lambda1 > 0) && (lambda2 > 0)
            equilibriumType = 'Unstable Equilibrium';
        else
            equilibriumType = 'Saddle Point';
        end
    
    else % Complex eigenvalues
        if real(lambda1) < 0 && real(lambda2) < 0
            equilibriumType = 'Stable Spiral';
        elseif real(lambda1) > 0 && real(lambda2) > 0
            equilibriumType = 'Unstable Spiral';
        else
            equilibriumType = 'Saddle Focus';
        end
    end
    
    % Print the result
    fprintf('The type of equilibrium is: %s \n', equilibriumType);
    fprintf('\n ------------------------------- \n');

end


% RATE OF FIRING  vs APPLIED CURRENT

rate_of_firing = zeros(21, 1);
applied_current = zeros(21, 1);
i = 1;

for Iext = 80:1:100
    
    syms V w
    Vnc1_eqn = (1/C)*(Iext - gCa*(0.5*(1+tanh((V-v1)/v2)))*(V-VCa) - gK*w*(V-VK) - gL*(V-VL)) == 0;
    wnc1_eqn = (0.5*(1+tanh((V-v3)/v4)) - w) == 0;
    eq_pt_1 = vpasolve([Vnc1_eqn, wnc1_eqn], [V, w]);

    V_eq_4 = double(eq_pt_1.V);
    w_eq_4 = double(eq_pt_1.w);

    tSpan = [0 1000];
    initial = [V_eq_4 + 0.2, w_eq_4 + 0.005];
    [t1, S1] = ode15s(@(t,S)morris_lecar(t,S), tSpan, initial);
    
    rate_of_firing(i) = firing_rate(t1, S1);
    applied_current(i) = Iext;
    i = i + 1;
end

figure;
plot(applied_current, rate_of_firing, 'Color', 'g', 'LineWidth', 1.5);
xlabel('Applied Current');
ylabel('Firing rate');
title('(9) Rate of firing Action Potentials VS Applied Current');
grid on;

% -------------------- DONE --------------------------------------
% ----------------------------------------------------------------

% ------------------- QUESTION 10 -------------------------------
% ---------------------------------------------------------------
fprintf("\n --------------------- QUESTION 10 ------------------ \n");

% New parameters for the MLE model

C = 20.0;
gCa = 4.0;
gK = 8.0;
gL = 2.0;
VCa = 120.0;
VK = -84;
VL = -60;
v1 = -1.2;
v2 = 18;
v3 = 12;
v4 = 17.4;
v5 = 12;
v6 = 7.4;
phi = 0.0667;
Iext = 30;

tSpan = [0 2000];

figure;
hold on
Vnc1 = @(V) (Iext - gCa .*(0.5 .*(1 + tanh((V-v1)./v2))) .* (V-VCa) - gL .*(V-VL))./(gK .* (V-VK));
wnc1 = @(V) (0.5 .* (1 + tanh((V-v3)./v4)));
fplot(@(V) Vnc1(V), [-80 100], 'Color', 'r');
fplot(@(V) wnc1(V), [-80 100], 'Color', 'b');
ylim([0 1]);
xlabel('V (mV)');
ylabel('w');

syms V w
Vnc_eqn = (1/C) .* (Iext - gCa .* (0.5 .* (1 + tanh((V - v1)/v2))) .* (V - VCa) - gK .* w .*(V - VK) - gL .* (V - VL)) == 0;
wnc_eqn = (0.5 .* (1 + tanh((V - v3)./v4)) - w) == 0;
V_guess = [-40.5, -20, 5];
w_guess = [0, 0.025, 0.275];
letter = ['A', 'B', 'C'];
V_arr = zeros(3, 1);
w_arr = zeros(3, 1);

for i = 1:3
    equi_pt = vpasolve([Vnc_eqn, wnc_eqn], [V_guess(i), w_guess(i)]);
    
    V_eq = double(equi_pt.V);
    w_eq = double(equi_pt.w);
    V_arr(i) = V_eq;
    w_arr(i) = w_eq;
    plot(V_eq, w_eq, 'ko', 'MarkerFaceColor', 'k');
    text(V_eq, w_eq, letter(i), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

    syms V w
    f1 = (1/C) .*(Iext - gCa .*(0.5 .*(1 + tanh((V - v1) ./ v2))) .* (V - VCa) - gK .* w .* (V - VK) - gL .* (V - VL));
    f2 = phi * ((0.5 .* (1 + tanh((V - v3) ./ v4))) - w) .* cosh((V - v3) ./ (2 .* v4));
    
    % JACOBIAN matrix
    J = jacobian([f1; f2], [V, w]);
    
    % Assuming V_eq and w_eq are already defined as equilibrium points
    % Calculate Jacobian at the equilibrium points
    J_equi = subs(J, {V, w}, {V_eq, w_eq});  
    J_equi = double(J_equi);
    
    % Calculate eigenvalues
    eigenvalues = eig(J_equi);
    
    lambda1 = eigenvalues(1);
    lambda2 = eigenvalues(2);
    
    % Print Jacobian and Eigenvalues    
    fprintf("Jacobian @ [%.4f, %.4f]: \n", V_eq, w_eq);  % Print equilibrium points
    disp(J_equi);
    disp("Eigenvalues: ")
    disp(eigenvalues);
    
    % Check the TYPE OF EQUILIBRIUM
    if isreal(lambda1) && isreal(lambda2) % Both eigenvalues are real
        if (lambda1 < 0) && (lambda2 < 0)
            equilibriumType = 'Stable Equilibrium';
        elseif (lambda1 > 0) && (lambda2 > 0)
            equilibriumType = 'Unstable Equilibrium';
        else
            equilibriumType = 'Saddle Point';
        end
    
    else % Complex eigenvalues
        if real(lambda1) < 0 && real(lambda2) < 0
            equilibriumType = 'Stable Spiral';
        elseif real(lambda1) > 0 && real(lambda2) > 0
            equilibriumType = 'Unstable Spiral';
        else
            equilibriumType = 'Saddle Focus';
        end
    end

    if strcmp(equilibriumType, 'Saddle Point')
        [eigVectors, D] = eig(J_equi);
    end
    
    % Print the result
    fprintf('Equilibrium for point %s is: %s\n', letter(i), equilibriumType);
    fprintf('\n ------------------------------------------- \n')
end
legend('V nullcline', 'w nullcline', 'equi point A', 'equi point B', 'equi point C');
title('(10) Phase Plane Plot for Iext = 30 [3 equilibrium points]');
grid on;
hold off;


figure;
hold on
Vnc1 = @(V) (Iext - gCa .*(0.5 .*(1 + tanh((V-v1)./v2))) .* (V-VCa) - gL .*(V-VL))./(gK .* (V-VK));
wnc1 = @(V) (0.5 .* (1 + tanh((V-v3)./v4)));
fplot(@(V) Vnc1(V), [-80 100], 'Color', 'r');
fplot(@(V) wnc1(V), [-80 100], 'Color', 'b');
for i = 1:3
    plot(V_arr(i), w_arr(i), 'ko', 'MarkerFaceColor', 'k');
end
ylim([0 1]);
xlabel('V (mV)');
ylabel('w');

% Plot the manifolds of the SADDLE POINT:
vs = eigVectors(:,1)';
vu = eigVectors(:,2)';
V_eq2 = V_arr(2);
w_eq2 = w_arr(2); 

% Evaluating manifolds at different initial points w.r.t to either side of
% the saddle node equilibrium point

start1 = [V_eq2, w_eq2] +  1 * vs;
start2 = [V_eq2, w_eq2] +  1 * vu;
start3 = [V_eq2, w_eq2] -  1 * vs;
start4 = [V_eq2, w_eq2] -  1 * vu;

[t1,S1]= ode15s(@(t,S)morris_lecar(t,S), [0 100], start1);
plot(S1(:,1), S1(:,2), 'm', 'LineWidth', 1.5, 'LineStyle', '--');

[t1,S1]= ode15s(@(t,S)morris_lecar(t,S), [0 -80], start2);
plot(S1(:,1), S1(:,2), 'g', 'LineWidth', 1.5, 'LineStyle', '--');

[t1,S1]= ode15s(@(t,S)morris_lecar(t,S), [0 100], start3);
plot(S1(:,1), S1(:,2), 'm', 'LineWidth', 1.5, 'LineStyle', '--');

[t1,S1]= ode15s(@(t,S)morris_lecar(t,S), [0 -30], start4);
plot(S1(:,1), S1(:,2), 'g', 'LineWidth', 1.5, 'LineStyle', '--');

grid on;
title('(10) Manifolds of the Saddle Node');

% ------------- DONE ----------------------------------------------

% ---------------------- QUESTION 11 ---------------------------------
% --------------------------------------------------------------------

% Nature of Equilibrium points for current range between 30 and 50 μA/cm2
% specially between 39 - 40 μA/cm2
gCa = 4;
VCa = 120;
gK = 8;
VK = -84;
gL = 2;
VL = -60;
v1 = -1.2;
v2 = 18;
v3 = 12;
v4 = 17.4;
phi = 0.0667;
C = 20;

currents = [30, 35, 39, 39.1, 39.2, 39.3, 39.4, 39.5, 39.5, 39.6 39.7, 39.8, 39.9, 40, 42.5, 45, 47.5, 50];
N = size(currents);
N = N(2);
frequency = zeros(N, 1);

fprintf("\n ------------- QUESTION 11 ------------------------ \n")
for i = 1:N
    % Plotting the Phase Plane with V and w null clines
    Iext = currents(i);
    
    if (Iext > 35) && (Iext < 41)       
        figure;
        hold on
        
        Vnc1 = @(V) (Iext - gCa .*(0.5 .*(1 + tanh((V-v1)./v2))) .* (V-VCa) - gL .*(V-VL))./(gK .*(V - VK));
        wnc1 = @(V) (0.5 .* (1 + tanh((V - v3)./v4)));
        fplot(@(V) Vnc1(V), [-80 100], '--');
        fplot(@(V) wnc1(V), [-80 100], '-');
        xlabel('V (in mV)');
        ylabel('w');
        grid on;

        plot_title = strcat('Phase Plane Plot of MLE for Iext   ', num2str(Iext));
        title(plot_title);
        saveas(gcf, [plot_title, '.jpg']);
        hold on;
    end
    
    % Nature of Equilibrium points
    syms V w
    Vnc1_eqn = (1/C) .* (Iext - gCa .* (0.5 .* (1 + tanh((V - v1) ./ v2))) .* (V - VCa) - gK .* w .* (V - VK) - gL .* (V - VL)) == 0;
    wnc1_eqn = (0.5 .* (1 + tanh((V-v3)./v4)) - w) == 0;
    soln1 = vpasolve([Vnc1_eqn, wnc1_eqn], [V, w], [-40, 0]);
    soln2 = vpasolve([Vnc1_eqn, wnc1_eqn], [V, w], [-20, 0]);
    soln3 = vpasolve([Vnc1_eqn, wnc1_eqn], [V, w], [0, 0.2]);

    V_eq1 = double(soln1.V);
    w_eq1 = double(soln1.w);
    
    V_eq2 = double(soln2.V);
    w_eq2 = double(soln2.w);
    
    V_eq3 = double(soln3.V);
    w_eq3 = double(soln3.w);
    
    % THERE CAN BE AT MOST 3 EQUILIBRIUM POINTS

    if isequal([V_eq1, w_eq1], [V_eq2, w_eq2]) 
        fprintf("For Iext = %f, One distinct equilibrium point \n", Iext);

        [t,S] = ode15s(@(t,S)morris_lecar(t,S), [0 3000], [V_eq3 + 0.5, w_eq3 + 0.05]);
        
        if (Iext > 35) && (Iext < 41)
            plot(V_eq3, w_eq3, 'bo');
            text(V_eq3, w_eq3, ['(' num2str(round(V_eq3,3)) ',' num2str(round(w_eq3,3)) ')']);
            grid on;
            plot(S(:,1), S(:,2));
        end
        
        frequency(i) = firing_rate(t, S);

    else

        fprintf("For Iext = %f, Three distinct equilibrium points \n", Iext);
        [t,S] = ode15s(@(t,S)morris_lecar(t,S), [0 3000], [V_eq2 - 0.1, w_eq2 - 0.01]);
        frequency(i) = 0;
        
        if (Iext >= 35) && (Iext < 41)
            plot(V_eq1, w_eq1, 'o');
            grid on;
            plot(V_eq2, w_eq2, '-k');
            grid on;
            plot(V_eq3, w_eq3, 'o');
            text(V_eq3, w_eq3, ['(' num2str(round(V_eq3,3)) ',' num2str(round(w_eq3,3)) ')']);
            grid on;
            plot(S(:,1), S(:,2));
        end
    end
end

figure;
hold on
ylabel('Firing Rate (in hz or 1/s)');
xlabel('Iext (in uA)');
grid on;
plot_title = strcat('Firing rate for Iext ', num2str(Iext));

title(plot_title);
saveas(gcf, [plot_title, '.jpg']);
plot(currents, frequency);
hold off    


% --------------- DONE --------------------------------------


% ------------- HODGKIN HUXLEY EQUATIONS ------------------------
% --------------------------------------------------------------

% -------------------------- QUESTION 12 --------------------------

function dS = hodgkin_huxley(t, S)
    global C; % Membrane capacitance
    global Iext; % Externally applied current

    global gK; % Potassium channel
    global gNa; % Sodium channel
    global gL; % Leak channel
    
    global VK;
    global VNa;
    global VL;
    
    global phi;
    global tol; % Tolerance term to avoid any 0/0 concitions
    
    V = S(1);
    n = S(2);
    m = S(3);
    h = S(4);
    tol = 1e-6;

    alpha_n =  - (0.01 .* phi .* (V + tol + 50)) ./ (exp(-(V + tol + 50)./10) - 1);
    alpha_m =  - (0.1 .* phi .* (V + tol + 35)) ./ (exp(-(V + tol + 35)./10) - 1);
    alpha_h = 0.07 .* phi .* exp(-(V + 60)./20);
    
    beta_n = 0.125 .* phi .* exp(-(V + 60)./80);
    beta_m = 4 .* phi .* exp(-(V + 60)./18);
    beta_h = phi ./ (exp(-(V + 30)./10) + 1);

    dV_dt = (1/C) .* (Iext - gK .* n^4 .* (V - VK) - gNa .* m^3 .* h .* (V - VNa) - gL .* (V - VL));
    dn_dt = alpha_n .* (1 - n) - beta_n .* n;
    dm_dt = alpha_m .* (1 - m) - beta_m .* m;
    dh_dt = alpha_h .* (1 - h) - beta_h .* h;
    
    dS = [dV_dt; dn_dt; dm_dt; dh_dt];
end

% ----------------------- QUESTION 13 ------------------------------
% ------------------------------------------------------------------
% Find the value of EL necessary to make the resting potential at -60 mV

global C gK VK gNa VNa gL Iext tol VL

tol = 1e-6;
Vr = - 60;
Iext = 0;
gK = 36.0; % mS/cm^2
gNa = 120.0; % mS/cm^2
gL = 0.3; % mS/cm^2
VK = - 72.0; % mV
VNa = 55.0; % mV
phi = 1; % Temperature coefficient

alpha_n =  - (0.01 .* phi .* (Vr + tol + 50)) ./ (exp(-(Vr + tol + 50)./10) - 1);
alpha_m =  - (0.1 .* phi .* (Vr + tol + 35)) ./ (exp(-(Vr + tol + 35)./10) - 1);
alpha_h = 0.07 .* phi .* exp(-(Vr + 60)./20);

beta_n = 0.125 .* phi .* exp(-(Vr + 60)./80);
beta_m = 4 .* phi .* exp(-(Vr + 60)./18);
beta_h = phi ./ (exp(-(Vr + 30)./10) + 1);

n_inf = alpha_n ./ (alpha_n + beta_n);
m_inf = alpha_m ./ (alpha_m + beta_m);
h_inf = alpha_h ./ (alpha_h + beta_h);

E_leak = Vr - (1/gL) .* (Iext - gK .* (n_inf^4) .* (Vr - VK) - gNa .* (m_inf)^3 .* h_inf .* (Vr - VNa));
fprintf('\n ------------------------- QUESTION 13 ------------------------------ \n');
fprintf('EL = %.3f mV \n', E_leak);
VL = E_leak;

% ------------------------------------------------------------ %
% ------------------------------------------------------------ %

C = 1.0;
gK = 36.0;
VK = -72.0;
gNa = 120.0;
VNa = 55.0;
gL = 0.3;
Iext = 0;
tol = 1e-10;

Iext = 10;
tSpan = [0 100];
initial = [-60, n_inf, m_inf, h_inf]; 
[t1, S] = ode15s(@(t,S)hodgkin_huxley(t,S), tSpan, initial, options);
figure;
plot(t1, S(:,1));
xlabel('Time (ms)');
ylabel('Membrane Voltage (mV)');
grid on;
title('Simulation of HH model with given parameters');    

% -----------------------------------------------------------------------

% ---------------------- QUESTION 14 --------------------------
% -------------------------------------------------------------

% STABILITY OF HH MODEL AT REST WITH Iext = 0
function equi_pt = get_equilibrium_point(Iext)
    global C;
    global gK;
    global gNa;
    global gL;
    global VK;
    global VNa;
    global VL;
    global tol;
    
    syms V n m h 
    alpha_n =  -0.01 * (V + tol + 50)/(exp(-(V + tol + 50)/10)-1);
    alpha_m =  -0.1 * (V + tol + 35)/(exp(-(V + tol + 35)/10)-1);
    alpha_h = 0.07 * exp(-(V + 60)/20);
    
    beta_n = 0.125 * exp(-(V + 60)/80);
    beta_m = 4 * exp(-(V + 60)/18);
    beta_h = 1/(exp(-(V + 30)/10) + 1);
    
    n_inf = alpha_n/(alpha_n + beta_n);
    m_inf = alpha_m/(alpha_m + beta_m);
    h_inf = alpha_h/(alpha_h + beta_h);
    
    % Equation of nullclines
    V_nc = (1/C) .* (Iext - gK .* n^4 .* (V - VK) - gNa .* m^3 .* h .* (V - VNa) - gL .* (V - VL)) == 0;
    n_nc = alpha_n .* (1-n) - beta_n .* n == 0 ;
    m_nc = alpha_m .* (1 - m) - beta_m .* m == 0 ;
    h_nc = alpha_h .* (1 - h) - beta_h .* h == 0 ;
    
    equi_pt = vpasolve([V_nc, n_nc, m_nc, h_nc], [V, n, m ,h]);
end

equi_pt = get_equilibrium_point(0);
V_eq1 = double(equi_pt.V);
n_eq1 = double(equi_pt.n);
m_eq1 = double(equi_pt.m);
h_eq1 = double(equi_pt.h);
  
fprintf("Equilibrium points for Iext = 0 \n");
fprintf("V = %.2f \n", V_eq1);
fprintf("n = %.2f \n", n_eq1);
fprintf("m = %.2f \n", m_eq1);
fprintf("h = %.2f \n", h_eq1);


% Threshold of the model for brief current pulses
Iext = 0;
impulse_arr = linspace(0, 20, 100);
max_V = zeros(100, 1);

for i = 1:100
    V_initial = -60 + impulse_arr(i)/C;
    initial = [V_initial, n_inf, m_inf, h_inf];
    [t, S] = ode15s(@(t,S)hodgkin_huxley(t,S), tSpan, initial, options);
    max_V(i) = max(S(:,1));
end

figure;
plot(impulse_arr, max_V);
hold on;
xlabel('Current impulse (μA/cm^2)');
ylabel('Maximum membrane voltage (mV)')
title('Max voltage VS current pulse');
grid on;

min = min(max_V(:));
max = max(max_V(:));
threshold_voltage = 0.5 * (min + max);

for i = 1:100
    if (max_V(i) >= threshold_voltage)
        threshold_current = impulse_arr(i);
        break;
    end
end

fprintf("----------------- QUESTION 14 ----------------------------")
fprintf("\n Threshold Current = %.2f μA/cm^2 \n", threshold_current);

% ----------------------- DONE ----------------------------------

% ------------------------- QUESTION 15 -----------------------------
% Stability analysis of equilibrium points for steady current injections

function stability_analysis(Iext)
    global C;
    global gK;
    global gNa;
    global gL;
    global VK;
    global VNa;
    global VL;
    global tol;
    
    tol = 1e-10;
    syms V n m h 
    alpha_n =  -0.01 * (V + tol + 50)/(exp(-(V + tol + 50)/10)-1);
    alpha_m =  -0.1 * (V + tol + 35)/(exp(-(V + tol + 35)/10)-1);
    alpha_h = 0.07 * exp(-(V + 60)/20);
    
    beta_n = 0.125 * exp(-(V + 60)/80);
    beta_m = 4 * exp(-(V + 60)/18);
    beta_h = 1/(exp(-(V + 30)/10) + 1);
    
    n_inf = alpha_n/(alpha_n + beta_n);
    m_inf = alpha_m/(alpha_m + beta_m);
    h_inf = alpha_h/(alpha_h + beta_h);
    
    % Equation of nullclines
    V_nc = (1/C) .* (Iext - gK .* n^4 .* (V - VK) - gNa .* m^3 .* h .* (V - VNa) - gL .* (V - VL)) == 0;
    n_nc = alpha_n .* (1-n) - beta_n .* n == 0 ;
    m_nc = alpha_m .* (1 - m) - beta_m .* m == 0 ;
    h_nc = alpha_h .* (1 - h) - beta_h .* h == 0 ;
    
    equi_pt = vpasolve([V_nc, n_nc, m_nc, h_nc], [V, n, m ,h]);
    
    V_eq1 = double(equi_pt.V);
    n_eq1 = double(equi_pt.n);
    m_eq1 = double(equi_pt.m);
    h_eq1 = double(equi_pt.h);
      
    fprintf("Equilibrium points: \n");
    fprintf("V = %.2f \n", V_eq1);
    fprintf("n = %.2f \n", n_eq1);
    fprintf("m = %.2f \n", m_eq1);
    fprintf("h = %.2f \n", h_eq1);

    dV_dt = (1/C) * (Iext - gK * n^4 * (V - VK) - gNa * m^3 * h* (V - VNa) - gL * (V - VL));
    dn_dt = alpha_n * (1-n) - beta_n * n;
    dm_dt = alpha_m * (1 - m) - beta_m * m ;
    dh_dt = alpha_h * (1 - h) - beta_h * h ;
    
    Jacobian_arr = jacobian([dV_dt, dn_dt, dm_dt, dh_dt],[V,n,m,h]);
    
    V = V_eq1;
    n = n_eq1;
    m = m_eq1;
    h = h_eq1;
    
    Jmatrix = zeros(4,4);
    
    Jmatrix(1,1) = subs(Jacobian_arr(1,1));
    Jmatrix(1,2) = subs(Jacobian_arr(1,2));
    Jmatrix(1,3) = subs(Jacobian_arr(1,3));
    Jmatrix(1,4) = subs(Jacobian_arr(1,4));
    Jmatrix(2,1) = subs(Jacobian_arr(2,1));
    Jmatrix(2,2) = subs(Jacobian_arr(2,2));
    Jmatrix(2,3) = subs(Jacobian_arr(2,3));
    Jmatrix(2,4) = subs(Jacobian_arr(2,4));
    Jmatrix(3,1) = subs(Jacobian_arr(3,1));
    Jmatrix(3,2) = subs(Jacobian_arr(3,2));
    Jmatrix(3,3) = subs(Jacobian_arr(3,3));
    Jmatrix(3,4) = subs(Jacobian_arr(3,4));
    Jmatrix(4,1) = subs(Jacobian_arr(4,1));
    Jmatrix(4,2) = subs(Jacobian_arr(4,2));
    Jmatrix(4,3) = subs(Jacobian_arr(4,3));
    Jmatrix(4,4) = subs(Jacobian_arr(4,4));

    eigenValues = eig(Jmatrix);

    fprintf("\n Eigenvalues: \n");
    
    for k = 1:length(eigenValues)
        if imag(eigenValues(k)) == 0
            % Display only the real part
            fprintf(' %.4f\n', real(eigenValues(k)));
        else
            % Display both real and imaginary parts if there is an imaginary component
            fprintf('  %.4f %+.4fi\n', real(eigenValues(k)), imag(eigenValues(k)));
        end
    end

    fprintf("Type of equilibrium: ");
    
    if all(real(eigenValues) < 0)
        disp('Stable Equilibrium');
    elseif all(real(eigenValues) > 0)
        disp('Unstable Equilibrium');
    elseif any(real(eigenValues) < 0) && any(real(eigenValues) > 0)
        disp('Saddle point.');
    elseif all(real(eigenValues) == 0)
        disp('Equilibrium is a center or requires further analysis.');
    end
end

% ---------------------------------------------------------------------
current_injection = 8:0.5:12;
fprintf("\n --------------- QUESTION 15 ------------------------ \n");
for i = 1:numel(current_injection)
    Iext = current_injection(i);
    fprintf("For Iext = %.2f \n", Iext); 
    stability_analysis(Iext);
    fprintf(" -------------------------------------------------------- \n");
end

% ----------------------------- DONE -----------------------------
% ----------------------------------------------------------------

% ------------------------- QUESTION 16 ----------------------------
% ------------------------------------------------------------------

function dS = HH_reduced__n(t, S, Iext)
    global C;
    global gK;
    global gNa;
    global gL;
    global VK;
    global VNa;
    global VL;
    global tol;

    C = 1.0; %  μF/cm^2
    gK = 36.0; % mS/cm^2
    gNa = 120.0; % mS/cm^2
    gL = 0.3; % mS/cm^2
    VK = - 72.0; % mV
    VNa = 55.0; % mV
    VL = -49; % mV
    phi = 1; % Temperature coefficient
    m_inf = 0.05; 
    h_inf = 0.60; 

    V = S(1);
    n = S(2);
    
    alpha_n =  - (0.01 .* phi .* (V + tol + 50)) ./ (exp(-(V + tol + 50)./10) - 1);
    beta_n = 0.125 .* phi .* exp(-(V + 60)./80);
    
    dV_dt = (1/C) .* (Iext - gK .* (n)^4 .* (V - VK) - gNa .* (m_inf)^3 .* (h_inf) .* (V - VNa) - gL .* (V - VL));
    dn_dt = alpha_n .* (1 - n) - beta_n * n;
    
    dS = [dV_dt; dn_dt];
end

current = 8:1:12;
colormap = ['r', 'g', 'b', 'm', 'k'];

for i = 1:numel(current)
    Iext = current(i);
    tSpan = [0 100];
    intial = [-60, m_inf]; 

    [t1, S] = ode15s(@(t,S)HH_reduced__n(t, S, Iext), tSpan, intial, options);
    
    figure;
    plot(t1, S(:,1), colormap(i));
    grid on;
    xlabel('Time (ms)');
    ylabel('Membrane Potential (mV)');

    title(sprintf("HH V-m Reduced model for Iext = %.2f", Iext));
end

% ------------------------------ DONE -----------------------------------


% ------------------ QUESTION 17 ----------------------------------------
% ------------------ ANODE BREAK EXCITATION with HH MODEL ---------------

% Applying a hyperpolarising current for 20 ms
tSpan1 = [0 20];
Iext = -3;
initial1 = [-60, n_inf, m_inf, h_inf];
[t1, S1] = ode15s(@(t,S)hodgkin_huxley(t,S), tSpan1, initial1, options);

n_after_Anode = S1(end,2);
h_after_Anode = S1(end,4);

% Find the final values of V, n, m, h
Iext = 0;
tSpan2 = [20 100];
initial2 = [S1(end,1), S1(end,2), S1(end,3), S1(end,4)];
[t2, S2] = ode15s(@(t,S)hodgkin_huxley(t,S), tSpan2, initial2, options);

total_time = [t1; t2];
total_V = [S1(:,1); S2(:,1)];
figure;
plot(total_time, total_V, 'r');
xlabel('time (ms)');
ylabel('Membrane potential (mV)');
title('Anode Break Excitation for HH model');
grid on;

% --------------------------- DONE ----------------------------------

% ------------------------ QUESTION 18 ---------------------------------
% ------------ ANODE BREAK EXCITATION in HH REDUCED MODEL --------------

function dS = HH_reduced__m(t, S, n_inf, h_inf, Iext)
    global C;
    global gK;
    global gNa;
    global gL;
    global VK;
    global VNa;
    global VL;
    global tol;

    C = 1.0; %  μF/cm^2
    gK = 36.0; % mS/cm^2
    gNa = 120.0; % mS/cm^2
    gL = 0.3; % mS/cm^2
    VK = - 72.0; % mV
    VNa = 55.0; % mV
    VL = -49; % mV
    phi = 1;
    tol = 1e-10;

    V = S(1);
    m = S(2);
    
    alpha_m =  - (0.1 .* phi .* (V + tol + 35)) ./ (exp(-(V + tol + 35)./10) - 1);
    beta_m = 4 .* phi .* exp(-(V + 60)./18);
    
    dV_dt = (1/C) .* (Iext - gK .* (n_inf)^4 .* (V - VK) - gNa .* (m)^3 .* (h_inf) .* (V - VNa) - gL .* (V - VL));
    dm_dt = alpha_m .* (1 - m) - beta_m * m;
    
    dS = [dV_dt; dm_dt];
end

function Characterize_equilibrium(V_eq, m_eq, n_Inf1, h_Inf1, Iext)
    syms V m

    C = 1.0; %  μF/cm^2
    gK = 36.0; % mS/cm^2
    gNa = 120.0; % mS/cm^2
    gL = 0.3; % mS/cm^2
    VK = - 72.0; % mV
    VNa = 55.0; % mV
    VL = -49; % mV
    tol = 1e-10;
    phi = 1; 

    alpha_m =  - (0.1 .* phi .* (V + tol + 35)) ./ (exp(-(V + tol + 35)./10) - 1);
    beta_m = 4 .* exp(-(V + 60)./18);
    
    f1 = (1/C) .* (Iext - gK .* n_Inf1^4 .* (V - VK) - gNa .* m^3 .* h_Inf1 .* (V - VNa) - gL .* (V - VL));
    f2 = alpha_m .* (1-m) - beta_m .* m;
    
    J = jacobian([f1; f2], [V, m]);
    
    % Assuming V_eq and w_eq are already defined as equilibrium points
    % Calculate Jacobian at the equilibrium points
    J_equi = subs(J, {V, m}, {V_eq, m_eq});  
    J_equi = double(J_equi);
    
    % Calculate eigenvalues
    eigenvalues = eig(J_equi);
    
    lambda1 = eigenvalues(1);
    lambda2 = eigenvalues(2);
    
    % Print Jacobian and Eigenvalues    
    fprintf("Jacobian @ [%.4f, %.4f]: \n", V_eq, m_eq);  
    % Print equilibrium points
    disp(J_equi);
    disp("Eigenvalues: ")
    disp(eigenvalues);
    
    % Check the TYPE OF EQUILIBRIUM
    if isreal(lambda1) && isreal(lambda2) % Both eigenvalues are real
        if (lambda1 < 0) && (lambda2 < 0)
            equilibriumType = 'Stable Equilibrium';
        elseif (lambda1 > 0) && (lambda2 > 0)
            equilibriumType = 'Unstable Equilibrium';
        else
            equilibriumType = 'Saddle Point';
        end
    
    else % Complex eigenvalues
        if real(lambda1) < 0 && real(lambda2) < 0
            equilibriumType = 'Stable Spiral';
        elseif real(lambda1) > 0 && real(lambda2) > 0
            equilibriumType = 'Unstable Spiral';
        else
            equilibriumType = 'Saddle Focus';
        end
    end

    fprintf("Type of equilibrium: %s \n", equilibriumType);
    fprintf("\n");
end

tol = 1e-6;
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
Vr = -60;

alpha_n = @(V) -0.01 .* (V + tol + 50)./(exp(-(V + tol + 50)./10)-1);
alpha_m = @(V) -0.1 .* (V + tol + 35)./(exp(-(V + tol + 35)./10)-1);
alpha_h = @(V) 0.07 .* exp(-(V + 60)./20);

beta_n = @(V) 0.125 .* exp(-(V + 60)./80);
beta_m = @(V) 4 .* exp(-(V + 60)./18);
beta_h = @(V) 1./(exp(-(V + 30)./10) + 1);

m_inf = alpha_m(Vr)./(alpha_m(Vr) + beta_m(Vr));
n_inf = alpha_n(Vr)./(alpha_n(Vr) + beta_n(Vr));
h_inf = alpha_h(Vr)./(alpha_h(Vr) + beta_h(Vr));

Iext = -3;
tSpan1 = [0 20];
initial1 = [-60, m_inf];
[t1, S1] = ode15s(@(t,S)HH_reduced__m(t,S,n_inf,h_inf,Iext), tSpan1, initial1, options);
Vnc1 = @(V) (((Iext - gK .* (V - VK) .* (n_inf)^4- gL .* (V - VL))./(gNa .* h_inf .* (V - VNa)))^(1/3));
mnc = @(V) alpha_m(V) ./ (alpha_m(V) + beta_m(V));

figure;
hold on;
grid on;
fplot(@(V) Vnc1(V), [-80 100], ':');
fplot(@(V) mnc(V), [-80 100], '--');

Iext = 0;
tSpan2 = [20 100];
initial2 = [S1(end,1), S1(end,2)];

n_Inf1 = n_after_Anode;
h_Inf1 = h_after_Anode;

[t2, S2] = ode15s(@(t,S)HH_reduced__m(t,S, n_Inf1, h_Inf1, Iext), tSpan2, initial2, options);
Vnc2 = @(V) ( (Iext - gK .* (V - VK) .* (n_inf)^4- gL .* (V-VL))./(gNa .* h_inf .* (V - VNa)))^(1/3);
fplot(@(V) Vnc2(V), [-80 100], '-*');

total_V = [S1(:,1); S2(:,1)];
total_m = [S1(:,2); S2(:,2)];                                                              

plot(total_V, total_m, 'r');
ylabel('m - Sodium gating variable');
xlabel('Membrane potential (mV)');
title('Phase Plane for anode break excitation');
grid on;

% Equilibrium points  
Iext = -3;
syms V m;

V_nc = (1/C) .* (Iext - gK .* n_inf^4 .* (V - VK) - gNa .* m^3 .* h_inf .* (V - VNa) - gL .* (V - VL)) == 0;
m_nc = alpha_m .* (1-m) - beta_m .* m == 0 ;

eq_pt = vpasolve([V_nc, m_nc], [V, m]);
V_eq1 = double(eq_pt.V);
m_eq1 = double(eq_pt.m);

fprintf("--------------------- QUESTION 18 -------------- \n");
fprintf("Equilibrium Point for n and h @ Resting potential:\n");
fprintf("V = %f \n", V_eq1);
fprintf("m = %f \n", m_eq1);

Characterize_equilibrium(V_eq1, m_eq1, n_inf, h_inf, Iext);

syms V m;

Iext = 0;
V_nc = (1/C) .* (Iext - gK .* (n_Inf1)^4 .* (V - VK) - gNa .* m^3 .* h_Inf1 .* (V - VNa) - gL .* (V - VL)) == 0;
m_nc = alpha_m .* (1-m) - beta_m .* m == 0 ;

eq_pt = vpasolve([V_nc, m_nc], [V, m]);
V_eq2 = double(eq_pt.V);
m_eq2 = double(eq_pt.m);

fprintf("\n");
fprintf("Equilibrium Point for n and h @ final values of anodal stimulus :\n");
fprintf("V = %f \n", V_eq2);
fprintf("m = %f \n", m_eq2);

Characterize_equilibrium(V_eq2, m_eq2, n_Inf1, h_Inf1, Iext);

% --------------------------- DONE ---------------------------------
% ------------------------------------------------------------------