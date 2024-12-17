clear; clc;

mu = [0.1, 1, 100, 1000]; % SET OF MU VALUES
time = [0 40]; % TIME SCALE FOR RUNNING THE SIMULATION

% SIMULATION PLOT & SIMULATION TIME COMPARISON

for i = 1:length(mu)
    dydt = @(t,y) [mu(i)*y(2); -y(1)./mu(i) + (1 - y(1).^2).*y(2)];
   
    tic;
    figure(1)
    [t, y] = ode45(dydt, time, [1; 0]); % Simulation time for ODE45
    time45 = toc;

    subplot(2, 2, i)
    plot(t, y(:,1), 'b', 'DisplayName', 'y');
    hold on;
    plot(t, y(:,2), 'Color', 'r', 'DisplayName', '1/\mu*dy/dt');
    title(['ODE45 soln. with \mu = ', num2str(mu(i))]);
    xlabel("time");
    ylabel("y(t) & y'(t)");
    legend('show');
    grid on;

    tic;
    figure(2)
    [t, y] = ode15s(dydt, time, [1; 0]); % Simulation time for ODE15s
    time15s = toc;

    subplot(2, 2, i)
    plot(t, y(:,1), 'b', 'DisplayName', 'y');
    hold on;
    plot(t, y(:,2), 'Color', 'r', 'DisplayName', '1/\mu*dy/dt');
    title(['ODE15s soln. with \mu = ', num2str(mu(i))]);
    xlabel("time");
    ylabel("y(t) & y'(t)");
    legend('show');
    grid on;

    if(time45 < time15s)
        fprintf("For mu = %.1f, Faster = ODE45 \n", mu(i));
    else
        fprintf("For mu = %.1f, Faster = ODE15s \n", mu(i));
    end

end


% PHASE PLANE PLOT DIAGRAM

figure(3)

for i = 1:length(mu)
    dydt = @(t,y) [y(2); mu(i) * (1 - y(1).^2) * y(2) - y(1)];
    [t, y] = ode45(dydt, time, [1; 0]);
    
    sgtitle('Phase Plane plots for different \mu values');

    subplot(2, 2, i)
    plot(y(:,1), y(:,2), 'b')
    title(['\mu = ', num2str(mu(i))]);
    xlabel("y");
    ylabel("y'/\mu");
    grid on;
end