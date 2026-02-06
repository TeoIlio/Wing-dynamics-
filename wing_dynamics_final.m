%% Dynamics of aeronautical structures/ Assignment 1 : Dynamics of aircraft wing 

% Assignment in the contect of "Dynamic of Aeronautical Structures" 
% Academic year : 2025-2026 
% Department of Mechanical Engineering and Aeronautics
% Theodoros Iliopoulos 1093547- Stylianos Apostolou 1093590

%% 
clc
clear 

%% Given system characteristics

h = 0.15 ; % blade height [m] (y-axis)
L = 8.9 ; % blade length [m] (x-axis)
wi = 2 ; % blade width [m] (z-zaxis)
E = 160*1e9 ; % elasticity modulus [Pa] 
v = 0.275 ; % poisson ratio [-]
rho = 1578 ; % density [kg/m^3] 

%% Calculation of needed characteristics 

A = wi*h ;  % Area 
G = E/(2*(1+v)) ;  % Shear Modulus  [Pa]
Iyy = wi*(h^3)/12 ; % second moment of area around y-axis [m^4] 
Izz = h*(wi^3)/12 ; % second moment of area around z-axis [m^4] 
Ipp = Iyy + Izz ; % mass moment of inertia [m^4]
c2 = -0.2*((wi/h)^-0.8784)+0.3386; % empirical correction coefficient 
J = c2 * (h^3) * wi ; % constant derived from system characteristics
 
%% System definition 

system.h = h ;
system.L = L ;
system.w = wi ;
system.E = E ;
system.v = v ; 
system.rho = rho ;
system.A = A ;
system.G = G ;
system.Iyy = Iyy ; 
system.Izz = Izz ; 
system.Ipp = Ipp ;
system.c2 = c2 ;
system.J = J ;

%% Finite Element Discretization
% 5 , 8 ,10 , 12 , 15 , 18 , 27 
%-----------------------
FE.N = 10 ; % number of elements
%----------------------
FE.eff_DOF = 3*FE.N ; % degrees of freedom with boundary conditions applied 
FE.DOF = 3*(FE.N+1) ; % DOF before the boundary conditions
FE.Le = system.L / FE.N ; %  element length (equal length for all FE)

%% Task 1 - Compute the natural frequencies and present the first 5 eigenmodes. Ensure convergence for the first 10 eigenfrequencies


[K,M,~] = get_system_matrix(system,FE) ; 
[Ke,Me] = get_element_matrix(system,FE);
[X,natural_freq,lambda] = solve_eigenvalue_problem(K,M) ; 
[w,bx,by]=get_plot(X,natural_freq,system,FE,1); 

 



%% Task 2 - Compute the equivalent dynamic system from the superposition of eigenmodes. (1st mode / first 3 modes). Compute the forced vibration at the free end of the wing.

modal_M = get_generalized_matrix(X,M) ; 
modal_K = get_generalized_matrix(X,K) ; 
%--------------------
% Selection of the number of modes for the superposition
%--------------------
modes_num = 3 ;
% --------------------

forcing_omega = 2*pi*(natural_freq(1)+natural_freq(2))/2 ; 
forcing_frequency = forcing_omega/(2*pi) ; 
%forcing_frequency = 20  ;  % for the 2nd case where the engine acts excitation
max_freq= max(forcing_frequency,natural_freq(modes_num));  % max frequency needed for timestep computation 
dt = 0.01/max_freq ;  % timestep 
total_time = 35*(1/max_freq) ; % amount of periods 
t = 0:dt:total_time ;

excitation = sin(forcing_omega*t) ;   % initiate only the sinusoidal function of the excitation, mangitude is found seperately
F0= zeros(FE.eff_DOF,1) ;  % magnitude 
% apply the force in the correct node 
if mod(FE.N,2) == 0 
    middle = FE.N/2 ;
else 
    middle = (FE.N/2) + 0.5 ;  
end 

force_acts_on = (middle)*3 -2 ;  % applied in the middle, works as an index 
F0(force_acts_on,1) = 0.1 ;   % magnitude at the proper node -DoF
modal_Fo = X' * F0 ;   % equivalent magnitude in the modal space 
modal_F = modal_Fo * excitation ; % excitation in the modal space 

% [modal_M]q_dot_dot + [modal_K]q = modal_F --> [I]q_dot_dot + [diag(lambda]q = Modal_F -reduced problem 
% u = X*q (physical space = X * modal space )
q_total_undamped = zeros(modes_num,length(modal_F(1,:))) ;  % rows = number of modes used , columns= discretized time 
q_free_undamped = zeros(modes_num,length(modal_F(1,:)))  ;  % solution of the homogeneous - complementrary solution 
q_forced_undamped = zeros(modes_num,length(modal_F(1,:)))  ;  % solution of the non-homogeneous - particular solution


for k = 1:modes_num  % find the time-varying functions- modal coefficient of each mode 
    
    [q_free_undamped(k,:),q_forced_undamped(k,:),q_total_undamped(k,:)] = solve_DE(t,modal_Fo(k),forcing_frequency,natural_freq(k),true,k) ; 
    
end 

u = X(:,1:modes_num) * q_total_undamped ;  % the displacement field in the physical space (functions of time)

% plot tip displacement 
figure ; 
plot(t(1:1:end),u(end-2,1:1:end),"LineStyle","-","Color","green","LineWidth",1)
grid on ;
xlabel("Time [seconds]")
ylabel("Displacement W [meters]")
title(sprintf("Response (W) for forced vibrarion (f = %.4f)",(forcing_frequency)))

% #### GIF - only for testing purposes/crosscheck  . ! Frames are
% notequivalent to real time so response shown is not in real time but
% governed by fps . !! It automatically saves it upon run time , dont run  it for long periods, it needs memory 
% close all 
% filename = 'giffy.gif';
% w_help = u(1:3:end-2,:) ; 
% w_actual = [ zeros(1,size(w_help,2)) ; w_help] ; 
% min_w = min(w_actual(:)) ; 
% max_w = max(w_actual(:)) ; 
% sf_gap = 0.15 * (max_w - min_w + eps);
% 
% x_points = 0:FE.Le:system.L; 
% x_points = x_points.' ; 
% figure; 
% set(gcf,'Position',[100 100 900 450]);
% 
% for i = 1:1:length(w_actual(1,:))
% 
%     plot(x_points,w_actual(:,i),'*','LineStyle','-','Color','red','LineWidth',2,'MarkerSize',8)
%     grid on ;
%     xlabel('X [m]')
%     ylabel(' W(x)  [m]')
%     title(sprintf("Wing deformation at t = %.4f seconds , forcing frequency = %.3f Hz",t(i),forcing_frequency))
% 
%     axis([min(x_points) max(x_points) (min_w-sf_gap) (max_w+sf_gap)]);  % fixed axes
%     drawnow;
%     % Capture and write frame to GIF
%     frame = getframe(gcf);
%     im = frame2im(frame);
%     [A,map] = rgb2ind(im,256);
% 
%     if i == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.04);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.04);
%     end
% end 
% 
% disp(['Saved GIF: ', filename]);

%% Task 03 - Damping Ce=0.05*[Ke] , compute response for the damped system 

alpha = 0.0005 ; 

damping_ratio = 2*pi*natural_freq * alpha / 2 ; 


q_total_damped = zeros(modes_num,length(modal_F(1,:))) ;  % rows = number of modes used , cols= discretized time 
q_free_damped = zeros(modes_num,length(modal_F(1,:)))  ;
q_forced_damped = zeros(modes_num,length(modal_F(1,:)))  ; 

for k = 1:modes_num
    
    [q_free_damped(k,:),q_forced_damped(k,:),q_total_damped(k,:)] = solve_DE_damped(t,modal_Fo(k),damping_ratio(k),forcing_frequency,natural_freq(k),false,k) ; 
    
end 

u_damped = X(:,1:modes_num) * q_total_damped ;  % the displacement field in the real space (functions of time)

% plots for the damped case 
figure ; 
plot(t(1:1:end),u_damped(end-2,1:1:end),"LineStyle","-","Color","k","LineWidth",1)
grid on ;
xlabel("Time [seconds]")
ylabel("Displacement W [meters]")
title(sprintf(" DAMPED Response (W) for forced vibrarion (f = %.4f)",(forcing_frequency)))

figure; 
plot(t(1:1:end),u(end-2,1:1:end),"LineStyle","-","Color","green","LineWidth",1)
hold on;
plot(t(1:1:end),u_damped(end-2,1:1:end),"LineStyle","-","Color","k","LineWidth",1)
grid on ;
ylabel("Displacement W")
xlabel("Time [seconds]")
title("Combined plot: Damped vs Undamped system's displacement (task 03)")
legend('Undamped','Damped') ;

%% Task 4 - Compute the poles of the system and the right-left eigenvalues, find the response of the damped system 
 
[~,~,C] = get_system_matrix(system,FE) ;


A_mat_1 = zeros(FE.eff_DOF) ;  % DOF x DOF
A_mat_2 = eye(FE.eff_DOF) ;
A_mat_3 = - M \ K ;   
A_mat_4 = - M \ C ;
A_mat = [ A_mat_1 , A_mat_2 ; A_mat_3 , A_mat_4] ;  % matrix used in the first order equivalent problem

% B_mat_1 = zeros(FE.eff_DOF,1) ;  not needed since forces are not taken into consideration 
% B_mat_2 = inv(M)  ; 
% B_mat = [B_mat_1 ; (B_mat_2 *F0) ]  ; 


% solve the right eigenvalue problem; Find right eigenvectors and eigenvalues
[X_right,eig_A] = eig(A_mat) ;  

lambda_T4 = diag(eig_A) ;   % eigenvalues as a single column , derived from the normal problem 

% solve the adjoint eigenvalue problem; Find the left eigenvectors and eigenvalues

[X_left,eig_A_adjoint] = eig(A_mat') ;  
lamnda_T4_adjoint = diag(eig_A_adjoint) ; 

yes = X_left.' * X_right;  % crosscheck that left and right eigevectors are orthonormal

% find the solution as a function of time 
x = FE.Le : FE.Le : system.L ; 
IC_trial = 0.0001*exp(x) ;
IC = zeros(2*FE.eff_DOF,1) ; % initial condition
IC(1:3:FE.eff_DOF) = IC_trial ;  


u_T4 = zeros(2*FE.eff_DOF,length(t)) ;  % init solution matrix 

for  i = 1:length(t) 
    u_T4(:,i) = X_right *(diag(exp(lambda_T4*t(i)))) * (X_left.' * IC) ;  % solution 
    
end 

% plot tip displacement
figure;
plot(t(1:1:end),u_T4(end-2-FE.eff_DOF,1:1:end),"LineStyle","-","Color","red","LineWidth",1)
grid on ;
ylabel("Displacement W")
xlabel("Time [seconds]")
title("Displacement for the damped case. Complex eigensolution (task 04)")

%%  task 05 and task 06  - Numerical integration 
tic;  % Start timer



dt1 = 1/(5*natural_freq(1)) ;  % the timestep is larger than the minimum limit , the solution breaks 
dt2 = 1/(5*natural_freq(end)) ; 
 
total_time_to_end = 50000/natural_freq(end) ;  % total sim time 
sim_time = 0:dt2:total_time_to_end ;   % discretzed time 
a1 = 1/(dt2^2) ; % auxilary 
a2 = 1/(2*dt2) ; % auxilary 
Fi = F0 * sin(forcing_omega*sim_time) ; % force 
modal_Fi = modal_Fo * sin(forcing_omega*sim_time) ;  % modal force 
modal_C = alpha * modal_K ; % modal damping matrix 

x_current = zeros(FE.eff_DOF,1) ;  % IC ,
x_previous = zeros(FE.eff_DOF,1) ; % IC 
x_current_m = zeros(FE.eff_DOF,1) ;  % IC for modal space 
x_previous_m = zeros(FE.eff_DOF,1) ; % IC for modal space

x_total = zeros(FE.eff_DOF,length(sim_time)) ; % init solution matrix 
x_modal_total = zeros(FE.eff_DOF,length(sim_time)) ; % init modal solution matrix

%Central differences 
for i=1:length(sim_time) 
  
 x_next =  (a1*M + a2*C)\( Fi(:,i) - (K - 2*a1*M)*x_current - (a1*M -a2*C)*x_previous  ) ;
 x_modal_next = (a1*modal_M + a2*modal_C)\( modal_Fi(:,i) - (modal_K - 2*a1*modal_M)*x_current_m - (a1*modal_M -a2*modal_C)*x_previous_m) ;
 x_total(:,i) = x_next  ; % update 
 x_modal_total(:,i) = x_modal_next ;  % update 

 if i ~= 1   % specia lcare for the 1st time step ,update
    x_previous = x_total(:,(i-1)) ;  
    x_previous_m = x_modal_total(:,(i-1)) ;
 else 
     x_previous = zeros(FE.eff_DOF,1) ; 
     x_previous_m = zeros(FE.eff_DOF,1) ;
 end
 x_current = x_total(:,i) ;
 x_current_m = x_modal_total(:,i) ;

end 

u_modal_total = X * x_modal_total ; % for modal space, results need to be multiplied with the modal matrix 

%plot tip displacement
figure;
plot(sim_time(1:1:end),x_total(FE.eff_DOF-2,1:1:end),"LineStyle","-","Color","red","LineWidth",1)
grid on ;
ylabel("Displacement W")
xlabel("Time [seconds]")
title("Displacement found via explicit integration (Task 05)")

 
figure;
plot(sim_time(1:1:end),u_modal_total(FE.eff_DOF-2,1:1:end),"LineStyle","-","Color","blue","LineWidth",1)
grid on ;
ylabel("Displacement W")
xlabel("Time [seconds]")
title("Displacement found via explicit integration with the use of generalized matrices (Task 06)")

figure;
plot(sim_time(1:1:end),u_modal_total(FE.eff_DOF-2,1:1:end),sim_time(1:1:end),x_total(FE.eff_DOF-2,1:1:end))
grid on ;
ylabel("Displacement W")
xlabel("Time [seconds]")
title("Displacement found via explicit integration ")
legend("Modal matrices", "Physical space matrices")
elapsedTime = toc;  % Stop timer and get elapsed time
fprintf('Elapsed time: %.6f seconds\n', elapsedTime);

%% task 07 - Aeroelasticity 

vel = 20;  % air velocity 
aoa_zero = 1;  % initial angle of attack 
a_n = 0.25 ;   % newmark a
d_n = 0.5 ;  % newmark delta 
 

dt_7 =  0.01; % time step 
sim = 0:dt_7:10;  % discretized time 

% init matrices 
R = zeros(FE.eff_DOF,length(sim)) ;   %  force 
u_7 = zeros(FE.eff_DOF,length(sim)) ;   % displacement 
u_dot_7 = zeros(FE.eff_DOF,length(sim)) ;  % velocity 
u_dot_dot_7 = zeros(FE.eff_DOF,length(sim)) ;  % acceleration

help1 = 1/(a_n*(dt_7^2)) ;  % auxilary 
help2 = 1/(a_n*(dt_7)) ;  % auxilary

R(:,1) = Air_load(vel,2,deg2rad(aoa_zero),u_7(:,1),3,(system.w*FE.Le),FE.N) ;  % initial aerodynamic loads 
u_dot_dot_7(:,1) = M \ ( R(:,1) - C*u_dot_7(:,1) - K*u_7(:,1) );  % initial acceleration 

% Newmark Method 
for i = 1:1:length(sim)-1
    
    LHS =  help1*M + d_n*help2*C + K  ;  % left hand side
    RHS = R(:,i) + M*(help1*u_7(:,i) + help2*u_dot_7(:,i) +(1/(2*a_n) - 1)*u_dot_dot_7(:,i)) + C*(d_n*help2*u_7(:,i) + ...
            ((d_n/a_n) -1)*u_dot_7(:,i) + ((d_n/a_n) -2)*(dt_7/2)*u_dot_dot_7(:,i) )  ;  % right hand side 

    u_7(:,i+1) = LHS \ RHS ; % next step displacements 
    u_dot_dot_7(:,i+1) = help1*(u_7(:,i+1)-u_7(:,i)) - help2*u_dot_7(:,i) - (1/(2*a_n) -1)*u_dot_dot_7(:,i) ;  % next step acceleration
    u_dot_7(:,i+1) = u_dot_7(:,i) + (1-d_n)*dt_7*u_dot_dot_7(:,i) +d_n*dt_7*u_dot_dot_7(:,i+1) ; % next step velocity

    R(:,i+1) = Air_load(vel,2,deg2rad(1),u_7(:,i+1),3,(system.w*FE.Le),FE.N) ;  % next step aerodynamic loads, calculated based on the displacements last found

end 

% plot 
figure('Position', [100, 100, 1200, 400]);

subplot(1,2,1)
plot(sim, u_7(end-2,:), 'LineStyle', '-', 'Color', 'red', 'LineWidth', 1.5)
xlabel('Time [s]')
ylabel('Vertical Displacement [m]')
title('BOX WING - Tip Vertical Displacement')
grid on;

subplot(1,2,2)
plot(sim, u_7(end-1,:), 'LineStyle', '-', 'Color', 'blue', 'LineWidth', 1.5)
xlabel('Time [s]')
ylabel('Rotation β_x [rad]')
title('BOX WING - Tip Rotation')
grid on;


% % ---### similar  GIF as in task 02 ####### 
% close all 
% filename = 'giffy.gif';
% w_help = u_7(1:3:end-2,:) ; 
% w_actual = [ zeros(1,size(w_help,2)) ; w_help] ; 
% min_w = min(w_actual(:)) ; 
% max_w = max(w_actual(:)) ; 
% sf_gap = 0.15 * (max_w - min_w + eps);
% 
% x_points = 0:FE.Le:system.L; 
% x_points = x_points.' ; 
% figure; 
% set(gcf,'Position',[100 100 900 450]);
% 
% for i = 1:1:length(w_actual(1,:))
% 
%     plot(x_points,w_actual(:,i),'*','LineStyle','-','Color','red','LineWidth',2,'MarkerSize',8)
%     grid on ;
%     xlabel('X [m]')
%     ylabel(' W(x)  [m]')
%     title(sprintf("Wing deformation at t = %.4f seconds , forcing frequency = %.3f Hz",t(i),forcing_frequency))
% 
%     axis([min(x_points) max(x_points) (min_w-sf_gap) (max_w+sf_gap)]);  % fixed axes
%     drawnow;
%     % Capture and write frame to GIF
%     frame = getframe(gcf);
%     im = frame2im(frame);
%     [A,map] = rgb2ind(im,256);
% 
%     if i == 1
%         imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.04);
%     else
%         imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.04);
%     end
% end 
% 
% disp(['Saved GIF: ', filename]);

%% Task 07 - I BEAM 

% I BEAM 

% Geometry characteristics
Bflange = 0.1*system.w;   % 10% of otiginal width 
Bweb    = 0.1*Bflange;  % 10% of the flange 
Hweb    = 0.7*system.h;  % 70% of total height 
Hflange = 0.15*system.h;   % 15% of total height (for each flange )

c2_Ibeam = -0.2*((Bflange/system.h)^-0.8784)+0.3386; % empirical correction coefficient 
%J_Ibeam = c2_Ibeam * (system.h^3) * Bflange ; % constant derived from system characteristics
J_Ibeam = (1/3)*(2*Bflange*Hflange^3 + Hweb*Bweb^3);   % constant 

d_Ibeam = (Hweb + Hflange)/2;        % distance of flange-web centroids

Iyy_Ibeam = ( Bflange*system.h^3 -(Bflange-Bweb)*(system.h -2*Hflange)^3)/12 ;   
Izz_Ibeam = ( 2*Hflange*Bflange^3 + Hweb*Bweb^3 )/12 ;    

A_Ibeam = 2*Bflange*Hflange + Hweb*Bweb ; % area 
Ipp_Ibeam = Iyy_Ibeam + Izz_Ibeam ;   % polar moment of inertia 

% init characteristics in a struct in order to properly input it in the functions requirring a struct input 
beam.Bflange = Bflange ;
beam.Hflange = Hflange  ;
beam.Bweb = Bweb  ;
beam.Hweb = Hweb ; 
beam.E = E ;
beam.v = v ; 
beam.rho = rho ;
beam.A = A_Ibeam ;
beam.G = G ;
beam.Iyy = Iyy_Ibeam ; 
beam.Izz = Izz_Ibeam ; 
beam.Ipp = Ipp_Ibeam ;
beam.c2 = c2_Ibeam ;
beam.J = J_Ibeam ;

% create the new stiffnessm,mass,damping matrices 
[K_I,M_I,C_I] = get_system_matrix(beam,FE) ; 
[Ke_I,Me_I] = get_element_matrix(beam,FE) ;

 
dt_7_I =  0.01;  % time step 
sim_I = 0:dt_7_I:10;  % discretized time step 

% init matrices 
R_I = zeros(FE.eff_DOF,length(sim_I)) ;
u_7_I = zeros(FE.eff_DOF,length(sim_I)) ; 
u_dot_7_I = zeros(FE.eff_DOF,length(sim_I)) ; 
u_dot_dot_7_I = zeros(FE.eff_DOF,length(sim_I)) ; 

help1 = 1/(a_n*(dt_7_I^2)) ; 
help2 = 1/(a_n*(dt_7_I)) ; 

R_I(:,1) = Air_load(vel,2,deg2rad(aoa_zero),u_7_I(:,1),3,(system.w*FE.Le),FE.N) ;  % initial aerodynamic loads 
u_dot_dot_7_I(:,1) = M \ ( R(:,1) - C*u_dot_7(:,1) - K*u_7(:,1) );  % initial acceleration 

% Newmarck method 
for i = 1:1:length(sim_I)-1

    
    LHS_I =  help1*M_I + d_n*help2*C_I + K_I  ;  % left hand side
    RHS_I = R_I(:,i) + M_I*(help1*u_7_I(:,i) + help2*u_dot_7_I(:,i) +(1/(2*a_n) - 1)*u_dot_dot_7_I(:,i)) + C_I*(d_n*help2*u_7_I(:,i) + ...
            ((d_n/a_n) -1)*u_dot_7_I(:,i) + ((d_n/a_n) -2)*(dt_7_I/2)*u_dot_dot_7_I(:,i) )  ;  % right hand side 

        % u_7(:,i+1) = ( help1*M + d_n*help2*C + K ) \ (R(:,i+1) + M*(help1*u_7(:,i) + help2*u_dot_7(:,i) +(1/(2*a_n) - 1)*u_dot_dot_7(:,i)) + C*(d_n*help2*u_7(:,i) + ...
        %     ((d_n/a_n) -1)*u_dot_7(:,i) + ((d_n/a_n) -2)*(dt_7/2)*u_dot_dot_7(:,i) ) );  

    u_7_I(:,i+1) = LHS_I \ RHS_I ; % next step displacements 
    u_dot_dot_7_I(:,i+1) = help1*(u_7_I(:,i+1)-u_7_I(:,i)) - help2*u_dot_7_I(:,i) - (1/(2*a_n) -1)*u_dot_dot_7_I(:,i) ;  % next step acceleration
    u_dot_7_I(:,i+1) = u_dot_7_I(:,i) + (1-d_n)*dt_7_I*u_dot_dot_7_I(:,i) +d_n*dt_7_I*u_dot_dot_7_I(:,i+1) ; % next step velocity

    R_I(:,i+1) = Air_load(vel,2,deg2rad(aoa_zero),u_7_I(:,i+1),3,(system.w*FE.Le),FE.N) ;  % next step aerodynamic loads, calculated based on the discplamenets last found

end 

figure('Position', [100, 100, 1200, 400]);

subplot(1,2,1)
plot(sim_I, u_7_I(end-2,:), 'LineStyle', '-', 'Color', 'red', 'LineWidth', 1.5)
xlabel('Time [s]')
ylabel('Vertical Displacement [m]')
title('I BEAM - Tip Vertical Displacement')
grid on;

subplot(1,2,2)
plot(sim_I, u_7_I(end-1,:), 'LineStyle', '-', 'Color', 'blue', 'LineWidth', 1.5)
xlabel('Time [s]')
ylabel('Rotation β_x [rad]')
title(' I BEAM - Tip Rotation')
grid on;



%% Element stiffness and mass matrices with Timoshenko assumptions
% Calculate the mass and stiffness element matrices based on the system's characteristics,as calculated/given initially 

% Inputs
%   -> system : A struct containing all characteristics needed to build the elemental matrices such as G , E, J, c2, A, rho etc 
%   -> FE : A struct containg all the information regarding the FE discretization including element number,length, dofs etc 
% Outputs:
%   -> Ke, Me : elemental matrices

function [Ke,Me] = get_element_matrix(system,FE)  

G = system.G ; 
A = system.A ;
Le = FE.Le ;
J = system.J ; 
E = system.E ; 
Iyy = system.Iyy ;
rho = system.rho ;
Ipp = system.Ipp; 

% element sitffness matrix 
Ke = [G*A/Le 0 -G*A/2 -G*A/Le 0 -G*A/2
    0 G*J/Le 0 0 -G*J/Le 0
    -G*A/2 0 (E*Iyy/Le + G*A*Le/3) G*A/2 0 (-E*Iyy/Le + G*A*Le/6)
    -G*A/Le 0 G*A/2 G*A/Le 0 G*A/2
    0 -G*J/Le 0 0 G*J/Le 0
    -G*A/2 0 (-E*Iyy/Le + G*A*Le/6) G*A/2 0 (E*Iyy/Le + G*A*Le/3)] ; 

% element mass matrix (conistent mass approach) 
Me = rho*[A*Le/3 0 0 A*Le/6 0 0
    0 Ipp*Le/3 0 0 Ipp*Le/6 0
    0 0 Iyy*Le/3 0 0 Iyy*Le/6
    A*Le/6 0 0 A*Le/3 0 0
    0 Ipp*Le/6 0 0 Ipp*Le/3 0
    0 0 Iyy*Le/6 0 0 Iyy*Le/3] ;

end

%% Build the system's Mass matrix and Stiffness matrix  
% each FE has 2 nodes , each node has 3 DOF(w,bx,by) apart from the 1st
% node that is clamped 

% Inputs : Same as the previous function 
% Output :  K, C, M matrices for the whole system 
% ** The C matrix is found for proportional damping assumption

function [K_total,M_total,C_total] = get_system_matrix(system,FE)


K_total = zeros(FE.DOF) ;
M_total = zeros(FE.DOF) ;
C_total = zeros(FE.DOF);

for element = 1:FE.N 
    i = element ;
    
    [Ke,Me] = get_element_matrix(system,FE) ;
    %map DOFs   
    map_start = 3*(i-1) + 1 ; % current element's 1st DOF position in the total matrices, effective DOF start position
    map_end = 3*(i) + 3 ;  % current element's 6th DOF position in the total matrices, effective DOF end position
    
    K_total(map_start:map_end,map_start:map_end) = K_total(map_start:map_end,map_start:map_end) + Ke ; 
    M_total(map_start:map_end,map_start:map_end) = M_total(map_start:map_end,map_start:map_end) + Me ; 
    C_total(map_start:map_end,map_start:map_end) = C_total(map_start:map_end,map_start:map_end) + 0.0005*Ke ; 
end
% exclude the 1st node (due to boundary conditions)
K_total(1:3,:) = [] ; 
K_total(:,1:3) = [] ; 

M_total(1:3,:) = [] ; 
M_total(:,1:3) = [] ; 

C_total(1:3,:) = [] ;
C_total(:,1:3) = [] ;

end 

%% find the eigenvalues,natural frequencies and eigenmodes

% Inputs 
function [X,natural_freq,lambda] = solve_eigenvalue_problem(K,M) 

%[K , M] = get_system_matrix(system,FE) ; 
[X,eig_D] = eig(K,M) ;  % x = eigenvector matrix(modal matrix) , eig_D = eigenvalues (diagonal matrix)

lambda = real(diag(eig_D)) ;
omega = sqrt(lambda) ; % rad /sec 
natural_freq = omega /(2*pi) ; % Hz 
[natural_freq,index] = sort(natural_freq,'ascend') ; %sort frequencies with ascending order 
X = X(:,index) ;   % modes are orthonormal to the M matrix 

end 

%% plot eigenmode 

% function made for visualizing modes 
% inputs:
%   -> X: the mode we want to visualize
%   -> natural_freq: the equivalent natural frequency 
%   -> system : struct, same as before
%   -> FE: struct , same as before
%   -> modes_nuumber :signal the mode that is being handled 
% Outputs:
%   -> Triple tile layout with modes for w,bx,by 

function  [w,bx,by]=get_plot(X,natural_freq,system,FE,modes_number)


w = X(1:3:FE.eff_DOF,:) ; 
bx = X(2:3:FE.eff_DOF,:) ; 
by = X(3:3:FE.eff_DOF,:) ;

x_points = 0:FE.Le:system.L; 
x_points = x_points.' ; 

for i=1:modes_number
    norm_w = w(:,i) / max(abs(w(:,i))) ; 
    norm_bx = bx(:,i) / max(abs(bx(:,i))) ; 
    norm_by = by(:,i) / max(abs(by(:,i))) ; 
    
    eff_w = [0;norm_w] ; 
    eff_bx = [0;norm_bx];
    eff_by = [0;norm_by];
    
    figure; 
    tiledlayout(3,1)
    % w(x) - x
    nexttile
    plot(x_points, eff_w, 'o--', 'Color','red', 'LineWidth',1)
    grid on ;
    xlabel('x [m]')
    ylabel('Displacement w(x)')
    title(sprintf('Mode: %d, frequency = %.4f Hz',i,natural_freq(i)))

    % bx(x) - x
    nexttile
    plot(x_points, eff_bx, 'o--', 'Color','blue', 'LineWidth',1)
    grid on;
    xlabel('x [m]')
    ylabel('Torsion bx(x)')
    title(sprintf('Mode: %d, frequency = %.4f Hz',i,natural_freq(i)))
    
    % by(x) - x
    nexttile
    plot(x_points, eff_by, 'o--', 'Color','green', 'LineWidth',1)
    grid on;
    xlabel('x [m]')
    ylabel('Bending by(x)')
    title(sprintf('Mode: %d, frequency = %.4f Hz',i,natural_freq(i)))

end 


end 

%% create the equivalent dynamic system in the modal space. Generalized mass /Generalized stiffness

% Generalized matrix creation 
% Inputs:
%   -> eigenmodes : The modal matrix 
%   -> matrix: the matrix for which we want to compute the modal equivalent
%   e.g. matrix = M -> finds the modal mass matrix 
% Outputs: 
%   -> generalized matrix 

function [modal_matrix] = get_generalized_matrix(eigenmodes,matrix)

modal_matrix = eigenmodes'* matrix * eigenmodes ; % generalized  matrix (modes are orthonomral to the mass matrix )
 % generalized stiffness matrix which equals to the eigenvalues diagonal matrix 

% modal_M and modal_K have, effectively, elements only in the diagonal. The off diagonal elements exist due to roundoff error. 
% generalized mass matrix is [I] and generalized stiffness is [lambda] diagonal

end

%% solve the 2nd order ode in the modal space  

% ODE Solver for the case of undamped reduced system (in modal space)
%    Inputs:
%       -> t: time vector 
%       -> modal_F_magntiude : the vector representing the force magnitude in modal space
%       -> forcing_frequency : the frequency of the excitation 
%       -> natural_frequency : the natural frequency for the mode we solve
%       -> show_plot : type Bool, true for showing plot, false for no plot
%       -> mode_flag : index for the mode being handled
%   Output:
%       -> complementary response (solution of homogeneous)
%       -> particular responser
%       -> total response ( addition of the above)

function [q_free,q_forced,q_total] = solve_DE(t,modal_F_magnitude,forcing_freqency,natural_frequency,show_plot,mode_flag)

%solve an undamed forced vibration with harmonic excitation 

% zero initial conditions 
q0 = 0  ; 
q0_dot = 0 ;
natural_omega = natural_frequency * 2*pi ;  % rad/sec
forcing_omega = forcing_freqency * 2*pi ;   % rad/sec

C = modal_F_magnitude/((natural_omega^2 - forcing_omega^2)) ; % coeff
A = q0 ;  % coeff
B = q0_dot/natural_frequency - C;  % coeff

q_free = A*cos(natural_omega*t) + B*sin(natural_omega*t) ; 
q_forced = C*sin(forcing_omega*t); 
q_total = q_free + q_forced ; 

% plots 
step = 1 ; 
if show_plot == true 
    figure;
    tiledlayout(3,1)
    
    nexttile
    plot(t(1:step:end),q_free(1:step:end),'--','Color','red','LineWidth',1)
    grid on ; 
    xlabel('Time [seconds]')
    ylabel(' Free Response ')
    title(sprintf('Free response -Mode: %d, natural frequency = %.4f Hz, forcing frequency = %.4f',mode_flag,natural_frequency,forcing_freqency))
    
    nexttile
    plot(t(1:step:end),q_forced(1:step:end),'--','Color','blue','LineWidth',1)
    grid on ; 
    xlabel('Time [seconds]')
    ylabel('Forced Response')
    title(sprintf(' Forced Response -Mode: %d, natural frequency = %.4f Hz ,forcing frequency = %.4f',mode_flag,natural_frequency,forcing_freqency))
    
    nexttile
    plot(t(1:step:end),q_total(1:step:end),'--','Color','green','LineWidth',1) 
    plot(t(1:step:end),q_total(1:step:end),'--','Color','green','LineWidth',1) 
    grid on ; 
    xlabel('Time [seconds]')
    ylabel('TotalResponse' )
    title(sprintf('Total Response- Mode: %d, natural frequency = %.4f Hz, forcing frequency = %.4f',mode_flag,natural_frequency,forcing_freqency))
end 

end 


%% solve 2nd order ode for damped system in modal space

% 2nd order ODE solver for damped system. Inputs and outputs are similar to
% the previous function 

function [q_free,q_forced,q_total] = solve_DE_damped(t,modal_F_magnitude,damping_ratio,forcing_freqency,natural_frequency,show_plot,mode_flag)

%solve an undamed forced vibration with harmonic excitation 

% zero initial conditions 
q0 = 0  ; 
q0_dot = 0 ;
natural_omega = natural_frequency * 2*pi ; 
forcing_omega = forcing_freqency * 2*pi ; 
damping_omega = natural_omega*sqrt(1-damping_ratio^2) ; 

helper = (forcing_omega/natural_omega)^2 ; % auxilary

C = (modal_F_magnitude / natural_omega^2) * ( (1- helper)/((1-helper)^2 +(2*damping_ratio*(sqrt(helper)))^2 ) ) ; 
D = (modal_F_magnitude / natural_omega^2) * ( (-2*damping_ratio*sqrt(helper))/((1-helper)^2 +(2*damping_ratio*(sqrt(helper)))^2 )) ;
A = -D ; 
B = (damping_ratio*natural_omega*A - C*forcing_omega )/damping_omega;

q_free = C*sin(forcing_omega*t) + D*cos(forcing_omega*t) ;  % particular
q_forced = exp(-damping_ratio*natural_omega*t) .* ( A*cos(damping_omega*t) + B*sin(damping_omega*t) ) ;  % complementary 
q_total = q_free + q_forced ; % actual response 

% plots 
step = 1 ; 
if show_plot == true 
    figure;
    tiledlayout(3,1)
    
    nexttile
    plot(t(1:step:end),q_free(1:step:end),'--','Color','red','LineWidth',1)
    grid on ; 
    xlabel('Time [seconds]')
    ylabel(' Free Response ')
    title(sprintf('(DAMPED)Free response -Mode: %d, natural frequency = %.4f Hz, forcing frequency = %.4f',mode_flag,natural_frequency,forcing_freqency))
    
    nexttile
    plot(t(1:step:end),q_forced(1:step:end),'--','Color','blue','LineWidth',1)
    grid on ; 
    xlabel('Time [seconds]')
    ylabel('Forced Response')
    title(sprintf(' (Damped)Forced Response -Mode: %d, natural frequency = %.4f Hz ,forcing frequency = %.4f',mode_flag,natural_frequency,forcing_freqency))
    
    nexttile
    plot(t(1:step:end),q_total(1:step:end),'--','Color','green','LineWidth',1) 
    grid on ; 
    xlabel('Time [seconds]')
    ylabel('TotalResponse' )
    title(sprintf('(Damped)Total Response- Mode: %d, natural frequency = %.4f Hz, forcing frequency = %.4f',mode_flag,natural_frequency,forcing_freqency))

end 

end 

