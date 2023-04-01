
% Define the start and end points of the trajectory
start_point1 = [0 0.5 0];
start_point = [0 0 1.8];
intermediate_point = [0.6 0 0.85];
end_point = [1.0 0 0.55];

% Define the time stamps for the trajectory
t = [0 1 2];

% Define the positions for the trajectory
p = [start_point; intermediate_point; end_point];

% Perform cubic interpolation
tq = linspace(0, 2, 4); % Define the query time stamps
pq = interp1(t, p, tq, 'cubic'); % Perform cubic interpolation

% Plot the interpolated trajectory
plot3(pq(:,1), pq(:,2), pq(:,3), 'LineWidth', 2);
xlabel('X'); ylabel('Y'); zlabel('Z');
grid on;
b=pq';