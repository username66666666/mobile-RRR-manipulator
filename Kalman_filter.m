%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%功能说明：Kalman在GPS中的应用
function Kalman_filter
T = 1;                               %雷达扫描周期
N = 71/T;                            %总的采样次数
X = zeros(4,N);                      %目标真实位置、速度
X(:,1) = [0,0,0,0];            %目标初始位置和速度
Z = zeros(4,N);                      %传感器对位置的量测
Z(:,1) = [X(1,1),X(1,2),X(1,3),X(1,4)];            %观测初始化
Q = diag([10,0.01,10,0.01]);     %过程噪声方差
R = [500,0,0,0;
    0,0.5,0,0;
    0,0,500,0;
    0,0,0,0.5];                      %观测噪声方差
A = [1,T,0,0;
    0,1,0,0;
    0,0,1,T;
    0,0,0,1]                        %一步转移阵
H = eye(4)                       %量测矩阵
I = eye(4);
B=[0.5*T^2, 0; T, 0;0,0.5*T^2;0,T];
u=[1,0.8];
u1=[1,0];
u2=[-1,-0.4]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 2:N-60
    X(:,k) = A*X(:,k-1) + B*u'+sqrtm(Q)*randn(4,1);   %目标真实轨迹
    Z(:,k)= H*X(:,k) + sqrtm(R)*randn(4,1);     %目标观测  
end

for k = N-59:N-50
    X(:,k) = A*X(:,k-1) + B*u1'+sqrtm(Q)*randn(4,1);   %目标真实轨迹
    Z(:,k)= H*X(:,k) + sqrtm(R)*randn(4,1);     %目标观测  
end

for k = N-49:N-20
    X(:,k) = A*X(:,k-1) +sqrtm(Q)*randn(4,1);   %目标真实轨迹
    Z(:,k)= H*X(:,k) + sqrtm(R)*randn(4,1);     %目标观测  
end

for k = N-19:N
    X(:,k) = A*X(:,k-1) +B*u2'+sqrtm(Q)*randn(4,1);   %目标真实轨迹
    Z(:,k)= H*X(:,k) + sqrtm(R)*randn(4,1);     %目标观测  
end
%Kalman滤波
X_kf = zeros(4,N);
X_kf(:,1) = X(:,1);                         %Kalman滤波状态初始化
P = eye(4);                                %协方差阵初始化
for k = 2:N-60
    X_pre = A*X_kf(:,k-1)+ B*u';                    %状态预测
    P_pre = A*P*A' + Q;              %协方差预测
    Kg = P_pre*H'*inv(H*P_pre*H' + R);        %计算Kalman增益
    e = Z(:,k) - H*X_pre;                     %新息
    X_kf(:,k) = X_pre + Kg*e;                 %状态更新
    P = (I-Kg*H)*P_pre;                %协方差更新
end

for k = N-59:N-50
    X_pre = A*X_kf(:,k-1)+ B*u1';                    %状态预测
    P_pre = A*P*A' + Q;              %协方差预测
    Kg = P_pre*H'*inv(H*P_pre*H' + R);        %计算Kalman增益
    e = Z(:,k) - H*X_pre;                     %新息
    X_kf(:,k) = X_pre + Kg*e;                 %状态更新
    P = (I-Kg*H)*P_pre;                %协方差更新
end

for k = N-49:N-20
    X_pre = A*X_kf(:,k-1);                    %状态预测
    P_pre = A*P*A' + Q;              %协方差预测
    Kg = P_pre*H'*inv(H*P_pre*H' + R);        %计算Kalman增益
    e = Z(:,k) - H*X_pre;                     %新息
    X_kf(:,k) = X_pre + Kg*e;                 %状态更新
    P = (I-Kg*H)*P_pre;                %协方差更新
end

for k = N-19:N
    X_pre = A*X_kf(:,k-1)+B*u2';                    %状态预测
    P_pre = A*P*A' + Q;              %协方差预测
    Kg = P_pre*H'*inv(H*P_pre*H' + R);        %计算Kalman增益
    e = Z(:,k) - H*X_pre;                     %新息
    X_kf(:,k) = X_pre + Kg*e;                 %状态更新
    P = (I-Kg*H)*P_pre;                %协方差更新
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure %画图显示
hold on;
box on;
plot(X(1,:),X(3,:),'-k');                    %真实轨迹
plot(Z(1,:),Z(3,:),'-b.');                   %观测轨迹
plot(X_kf(1,:),X_kf(3,:),'-r+');             %kalman滤波轨迹
legend('真实轨迹','观测轨迹','Kalman滤波轨迹');
xlabel('X/m');
ylabel('Y/m');
figure %画图显示
hold on;
box on;
plot(X(2,:),X(4,:),'-k'); %真实v
plot(Z(2,:),Z(4,:),'-b');
plot(X_kf(2,:),X_kf(4,:),'-r+');             %kalman滤波v
legend('真实v','measurement v','Kalman滤波轨迹');
xlabel('X:m/s');
ylabel('Y:m/s');

figure %画图显示 Y trajectory
hold on;
box on;
plot(0:N-1,X(3,:),'-k');                    %真实轨迹
plot(0:N-1,Z(3,:),'-b.');                   %观测轨迹
plot(0:N-1,X_kf(3,:),'-r+');             %kalman滤波轨迹
legend('真实轨迹','观测轨迹','Kalman滤波轨迹');
xlabel('t');
ylabel('Y/m');

figure %画图显示Y 速度
hold on;
box on;
plot(0:N-1,X(4,:),'-k'); 
plot(0:N-1,Z(4,:),'-b'); %真实轨迹
plot(0:N-1,X_kf(4,:),'-r+');             %kalman滤波轨迹
legend('真实轨迹','Kalman滤波轨迹');
xlabel('t');
ylabel('Y:m/s');

figure %画图显示 X trajectory
hold on;
box on;
plot(0:N-1,X(1,:),'-k');                    %真实轨迹
plot(0:N-1,Z(1,:),'-b.');                   %观测轨迹
plot(0:N-1,X_kf(1,:),'-r+');             %kalman滤波轨迹
legend('真实轨迹','观测轨迹','Kalman滤波轨迹');
xlabel('t');
ylabel('X/m');

figure %画图显示X velocity
hold on;
box on;
plot(0:N-1,X(2,:),'-k');%真实轨迹
plot(0:N-1,Z(2,:),'-b');
plot(0:N-1,X_kf(2,:),'-r+');             %kalman滤波轨迹
legend('真实轨迹','measurement v','Kalman滤波轨迹');
xlabel('t');
ylabel('X:m/s');
