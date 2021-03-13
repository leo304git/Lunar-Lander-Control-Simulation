%% Lander with AMEID and control - Numerical caculation
%Date:2020/11/19
%Code writter: Yi-Lun Hsu <r07921067@ntu.edu.tw> 
%Modified by Wei-Che Tseng <b07901165@ntu.edu.tw>
%Professor: Cheng-Wei Chen <cwchenee@ntu.edu.tw>
%All rights reserved.
%----------------------------------------------------------------------------------

clear;clc;
%% initial states setting
qi =[0 6 0.8 0.8 30*pi/180];        %I.C. of generalized coordinates: qi = [x(to) z(to) l1(to) l2(to) theta(to)];
ui =[0 -5 0 0 0];          %I.C. of generalized velocities:  ui = [x_dot(to) z_dot(to) l1_dot(to) l2_dot(to) theta_dot(to)];
xrmi  =[0 0 0];        %I.C. of AMEID states: xmi = [iar xpr xpr_dot]
xlmi  =[0 0 0];        %I.C. of AMEID states: xmi = [ial xpl xpl_dot]
lenq=length(qi);
%% PID parameters setting-----The digital P-controller is implemented using the so-called velocity form.
Kp= 500;
Kd = 12;
%Kp = 0;
%Kd = 0;
%% Parameters setting
%% 1-D Lander parameters
mb = 10;       %(kg)- Mass of the lander body 
m1 = 0.5;      %(kg)- Mass of leg
m2 = 0.5;
mp = 0.2;      %(kg)- Mass of VCM
mc = mp + mb;  %(kg)- Mass of mp+mb
Ig = 3.7309 %(kg/m^2)- Moment of Inertia of the lander
k = 7000;      %(N/m)- Buffer spring stiffness
c = k/20;      %(N.s/m)- Buffer damper coefficient
l1ss = 0.8;     %(m)- Length of the buffer
l2ss = 0.8;
W = 1;         %(m)- Width of the body
H = 1;         %(m)- Height of the body

S = 1.3856;
alpha = pi/3;
%% Gravity & Ground Properties
g = 9.80665;           %(m/s^2)-Standard gravity-the nominal gravitational acceleration of an object in a vacuum near the surface of the Earth
kf = 7.5e4;            %(N/m)-Ground spring stiffness
cf = 130;              %(N.s/m)-Ground damper coefficient
uk = 0.3;              %kinetic friction coefficient
%% AMEID parameter values
md=0.5;                %(kg)
mp=0.2;                %(kg)
La=6.4e-3;             %(H)
Ra=5.2;                %(Ohm)
kv=17;                 %(V.s/m)
kF=17;                 %(N/A)
L_stroke=0.5;          %(m)
%% Ameid state-space matrix
%pair(A,B,C)
Aam = [-Ra/La 0 -kv/La;0 0 1;kF/(mp+md) 0 0];
Bam = [1/La;0;0];
Cam = [kF 0 0];
%% Simulation condition
T = 5;                     %Simulation time 
dt = 0.0005;               %Sampling time
N = floor(T/dt);           %Steps
t = (0:1:N)'*dt;
H_Launch_delay = 700;
%% Definition of lander states
qn = zeros(N+1,lenq);      %Generalized coordinates
un = zeros(N+1,lenq);      %Generalized velocities

q1n = zeros(N+1,2);        %Footpad Mass coordinates
u1n = zeros(N+1,2);        %Footpad Mass velocities
q2n = zeros(N+1,2);        %Footpad Mass coordinates
u2n = zeros(N+1,2);        %Footpad Mass velocities

F1_impactn = zeros(N+1,1);  %Footpad1 impact force
F2_impactn = zeros(N+1,1);  %Footpad2 impact force
%% Definition of AMEID states
xmr = zeros(N+1,3);
xml = zeros(N+1,3);
%Input
Vr = zeros(N+1,1);
Vl = zeros(N+1,1);
%PID error matrix
ekr = zeros(N+2,1);
ekl = zeros(N+2,1);
%AMEID Force matrix
Frmeid = zeros(N+1,1);
Flmeid = zeros(N+1,1);

%setting initial condition
qn(1,:) = qi;              
un(1,:) = ui;              
xmr(1,:) = xrmi;
xml(1,:) = xlmi;
Launchr = false;
Launchl = false;

numr = 0;
numl = 0;
countr = 0;
countl = 0;
tic;

for i = 1:1:N
xbf = qn(i,1);
zbf = qn(i,2);
l1f = qn(i,3);
l2f = qn(i,4);
thetaf = qn(i,5);

xb_dotf = un(i,1);
zb_dotf = un(i,2);
l1_dotf = un(i,3);
l2_dotf = un(i,4);
theta_dotf = un(i,5);

q1n(i,:) = [ xbf - cos(thetaf)*(W/2 + S*cos(alpha)) - sin(thetaf)*(H/2 - S*sin(alpha)) - l1f*cos(alpha)*cos(thetaf) + l1f*sin(alpha)*sin(thetaf),...
                   zbf + cos(thetaf)*(H/2 - S*sin(alpha)) - sin(thetaf)*(W/2 + S*cos(alpha)) - l1f*cos(alpha)*sin(thetaf) - l1f*sin(alpha)*cos(thetaf)...
                   ];
q2n(i,:) = [ xbf + cos(thetaf)*(W/2 + S*cos(alpha)) - sin(thetaf)*(H/2 - S*sin(alpha)) + l2f*cos(alpha)*cos(thetaf) + l2f*sin(alpha)*sin(thetaf),...
                   zbf + cos(thetaf)*(H/2 - S*sin(alpha)) + sin(thetaf)*(W/2 + S*cos(alpha)) + l2f*cos(alpha)*sin(thetaf) - l2f*sin(alpha)*cos(thetaf)...
                   ];

u1n(i,:) = [ xb_dotf - l1_dotf*cos(alpha + thetaf) - theta_dotf*cos(thetaf)*(H/2 - S*sin(alpha)) + theta_dotf*sin(thetaf)*(W/2 + S*cos(alpha)) + theta_dotf*sin(alpha + thetaf)*l1f,...
                   zb_dotf - l1_dotf*sin(alpha + thetaf) - theta_dotf*cos(thetaf)*(W/2 + S*cos(alpha)) - theta_dotf*sin(thetaf)*(H/2 - S*sin(alpha)) - theta_dotf*cos(alpha + thetaf)*l1f...
                 ];
u2n(i,:) = [ xb_dotf + l2_dotf*cos(alpha - thetaf) + theta_dotf*sin(alpha - thetaf)*l2f - theta_dotf*cos(thetaf)*(H/2 - S*sin(alpha)) - theta_dotf*sin(thetaf)*(W/2 + S*cos(alpha)),...
                   zb_dotf - l2_dotf*sin(alpha - thetaf) + theta_dotf*cos(alpha - thetaf)*l2f + theta_dotf*cos(thetaf)*(W/2 + S*cos(alpha)) - theta_dotf*sin(thetaf)*(H/2 - S*sin(alpha))...
                 ];
    
gNb = qn(i,2);
gN1 = q1n(i,2);
gN2 = q2n(i,2);

rNb = un(i,2);
rN1 = u1n(i,2);  
rN2 = u2n(i,2);

gTb = qn(i,1);
gT1 = q1n(i,1);
gT2 = q2n(i,1);

rTb = un(i,1);
rT1 = u1n(i,1);
rT2 = u2n(i,1);

Mn = [m1 + m2 + mb,       0,   -m1*cos(alpha + thetaf), m2*cos(alpha - thetaf),... %1[1,2,3,4]
           (m1*(2*l1f*sin(alpha + thetaf) - 2*cos(thetaf)*(H/2 - S*sin(alpha)) + 2*sin(thetaf)*(W/2 + S*cos(alpha))))/2 - (m2*(2*cos(thetaf)*(H/2 - S*sin(alpha)) - 2*l2f*sin(alpha - thetaf) + 2*sin(thetaf)*(W/2 + S*cos(alpha))))/2;...
            0,      m1 + m2 + mb, -m1*sin(alpha + thetaf), -m2*sin(alpha - thetaf),... %2[1,2,3,4]
            (m2*(2*l2f*cos(alpha - thetaf) + 2*cos(thetaf)*(W/2 + S*cos(alpha)) - 2*sin(thetaf)*(H/2 - S*sin(alpha))))/2 - (m1*(2*l1f*cos(alpha + thetaf) + 2*cos(thetaf)*(W/2 + S*cos(alpha)) + 2*sin(thetaf)*(H/2 - S*sin(alpha))))/2;...
           -m1*cos(alpha + thetaf),  -m1*sin(alpha + thetaf), m1, 0,... %3[1,2,3, 4]
            (m1*(H*cos(alpha) + W*sin(alpha)))/2;... 
            m2*cos(alpha - thetaf), -m2*sin(alpha - thetaf), 0, m2,... %4[1,2,3,4]
            -(m2*(H*cos(alpha) + W*sin(alpha)))/2;... 
            m1*(l1f*sin(alpha + thetaf) - cos(thetaf)*(H/2 - S*sin(alpha)) + sin(thetaf)*(W/2 + S*cos(alpha))) - m2*(cos(thetaf)*(H/2 - S*sin(alpha)) - l2f*sin(alpha - thetaf) + sin(thetaf)*(W/2 + S*cos(alpha))),...
            m2*((W*cos(thetaf))/2 - (H*sin(thetaf))/2 + S*cos(alpha - thetaf) + l2f*cos(alpha - thetaf)) - m1*(S*cos(alpha + thetaf) + l1f*cos(alpha + thetaf) + (W*cos(thetaf))/2 + (H*sin(thetaf))/2),...
            (m1*(H*cos(alpha) + W*sin(alpha)))/2,...
            -(m2*(H*cos(alpha) + W*sin(alpha)))/2,...
            (H^2*m1)/4 + (H^2*m2)/4 + S^2*m1 + S^2*m2 + (W^2*m1)/4 + (W^2*m2)/4 + l1f^2*m1 + l2f^2*m2 + 2*S*l1f*m1 + 2*S*l2f*m2 + S*W*m1*cos(alpha) + S*W*m2*cos(alpha) - H*S*m1*sin(alpha) - H*S*m2*sin(alpha) + W*l1f*m1*cos(alpha) + W*l2f*m2*cos(alpha) - H*l1f*m1*sin(alpha) - H*l2f*m2*sin(alpha) + 37309/10000 ...
             ];
hn = [ 0;...
            -g*(m1 + m2 + mb);...
            g*m1*sin(alpha + thetaf) - (k*(2*l1f - 2*l1ss))/2 - c*l1_dotf;...
            g*m2*sin(alpha - thetaf) - (k*(2*l2f - 2*l2ss))/2 - c*l2_dotf;...
            g*m1*(S*cos(alpha + thetaf) + l1f*cos(alpha + thetaf) + (W*cos(thetaf))/2 + (H*sin(thetaf))/2) - g*m2*((W*cos(thetaf))/2 - (H*sin(thetaf))/2 + S*cos(alpha - thetaf) + l2f*cos(alpha - thetaf))...
          ];
%{
    hn = [                                       0;...
                                 - g*m1 - g*mc;...
      g*m1 - c*l1_dotf - (k*(2*l1f - 2*lss))/2];
 %}
LambdaN1 = kf*gN1*min(sign(gN1),0) - cf*rN1*min(sign(rN1),0)*min(sign(gN1),0);  
LambdaN2 = kf*gN2*min(sign(gN2),0) - cf*rN2*min(sign(rN2),0)*min(sign(gN2),0);
LambdaNb = kf*gNb*min(sign(gNb),0) - 10*cf*rNb*min(sign(rNb),0)*min(sign(gNb),0);  
LambdaN = [LambdaN1;LambdaN2;LambdaNb]; 
% 
LambdaT1 = -uk*LambdaN1*sign(rT1);
LambdaT2 = -uk*LambdaN2*sign(rT2);
LambdaTb = -uk*LambdaNb*sign(rTb);
LambdaT = [LambdaT1;LambdaT2;LambdaTb];
%
WN=[0, 0, 0;...
          1, 1, 1;...
          - cos(alpha)*sin(thetaf) - sin(alpha)*cos(thetaf), 0, 0;...
          0, cos(alpha)*sin(thetaf) - sin(alpha)*cos(thetaf), 0;...
          l1f*sin(alpha)*sin(thetaf) - sin(thetaf)*(H/2 - S*sin(alpha)) - l1f*cos(alpha)*cos(thetaf) - cos(thetaf)*(W/2 + S*cos(alpha)),...
          cos(thetaf)*(W/2 + S*cos(alpha)) - sin(thetaf)*(H/2 - S*sin(alpha)) + l2f*cos(alpha)*cos(thetaf) + l2f*sin(alpha)*sin(thetaf),...
          0];
WT=[1, 1, 1;...
         0, 0, 0;...
         sin(alpha)*sin(thetaf) - cos(alpha)*cos(thetaf), 0, 0;...
         0, sin(alpha)*sin(thetaf) + cos(alpha)*cos(thetaf), 0;...
         sin(thetaf)*(W/2 + S*cos(alpha)) - cos(thetaf)*(H/2 - S*sin(alpha)) + l1f*cos(alpha)*sin(thetaf) + l1f*sin(alpha)*cos(thetaf),...
         l2f*sin(alpha)*cos(thetaf) - sin(thetaf)*(W/2 + S*cos(alpha)) - l2f*cos(alpha)*sin(thetaf) - cos(thetaf)*(H/2 - S*sin(alpha)),...
         0];
P = WN*LambdaN+WT*LambdaT;
%%%% ------------------AMEID force------------------------%%%%
if Launchr == false
   if xmr(i,2)<0
      xmr(i,2)=0;
   end
   xmr(i+1,:) = ((Aam*dt+eye(3))*(xmr(i,:))')';
   hr = [ 0;...
             0;...
             0;...
             0;...
             0];
elseif Launchr == true
   countr = countr + 1;
   ekr(i+2) = qn(i,5);
   Vr(i+1) = Vr(i) + Kp*(ekr(i+2)-ekr(i+1)) + Kd*((ekr(i+2)-ekr(i+1))-(ekr(i+1)-ekr(i)))/dt;
   if xmr(i,2)<0
      xmr(i,2)=0;
   end
   xmr(i+1,:) = ((Aam*dt+eye(3))*(xmr(i,:))'+Bam*dt*Vr(i+1))';    
   hr = [ xmr(i,1)*kF*sin(thetaf);...
         -xmr(i,1)*kF*cos(thetaf);...
         0;...
         0;...
         -xmr(i,1)*kF*(W + cos(alpha)*(S + l1f))...
        ];
end

if Launchl == false
   if xml(i,2)<0
      xml(i,2)=0;
   end
   xml(i+1,:) = ((Aam*dt+eye(3))*(xml(i,:))')';
   hl = [ 0;...
             0;...
             0;...
             0;...
             0];
elseif Launchl == true
   countl = countl + 1;
   ekl(i+2) = -qn(i,5);
   Vl(i+1) = Vl(i) + Kp*(ekl(i+2)-ekl(i+1)) + Kd*((ekl(i+2)-ekl(i+1))-(ekl(i+1)-ekl(i)))/dt;;
   
   if xml(i,2)<0
      xml(i,2)=0;
   end
   xml(i+1,:) = ((Aam*dt+eye(3))*(xml(i,:))'+Bam*dt*Vl(i+1))';    
   hl = [ xml(i,1)*kF*sin(thetaf);...
         -xml(i,1)*kF*cos(thetaf);...
         0;...
         0;...
         xml(i,1)*kF*(W + cos(alpha)*(S + l2f))...
        ];
end
%-------L_stroke constraint------
if (xmr(i,2) >= L_stroke)
   xmr(i,2) = L_stroke;
   xmr(i+1,3) = 0;
   Launchr = false;
end
% || (numr >= 1)  numr = numr + 1;
if (xml(i,2) >= L_stroke)
   xml(i,2) = L_stroke;
   xml(i+1,3) = 0;
   Launchl = false;
end
 %|| (numl >= 1) numl = numl + 1;
%--------------------------------
%hu = [Fx;Fz;Fl];

%---Euler's method (Solve DAE)---
un(i+1,:) = un(i,:) + (Mn\(hn+P+hr+hl)*dt)';
qn(i+1,:) = qn(i,:) + un(i,:)*dt;
%--------------------------------

xbf = qn(i+1,1);
zbf = qn(i+1,2);
l1f = qn(i+1,3);
l2f = qn(i+1,4);
thetaf = qn(i+1,5);

xb_dotf = un(i+1,1);
zb_dotf = un(i+1,2);
l1_dotf = un(i+1,3);
l2_dotf = un(i+1,4);
theta_dotf = un(i+1,5);

u1n(i+1,:) = [ xb_dotf - l1_dotf*cos(alpha + thetaf) - theta_dotf*cos(thetaf)*(H/2 - S*sin(alpha)) + theta_dotf*sin(thetaf)*(W/2 + S*cos(alpha)) + theta_dotf*sin(alpha + thetaf)*l1f,...
                   zb_dotf - l1_dotf*sin(alpha + thetaf) - theta_dotf*cos(thetaf)*(W/2 + S*cos(alpha)) - theta_dotf*sin(thetaf)*(H/2 - S*sin(alpha)) - theta_dotf*cos(alpha + thetaf)*l1f...
                 ];
u2n(i+1,:) = [ xb_dotf + l2_dotf*cos(alpha - thetaf) + theta_dotf*sin(alpha - thetaf)*l2f - theta_dotf*cos(thetaf)*(H/2 - S*sin(alpha)) - theta_dotf*sin(thetaf)*(W/2 + S*cos(alpha)),...
                   zb_dotf - l2_dotf*sin(alpha - thetaf) + theta_dotf*cos(alpha - thetaf)*l2f + theta_dotf*cos(thetaf)*(W/2 + S*cos(alpha)) - theta_dotf*sin(thetaf)*(H/2 - S*sin(alpha))...
                 ];

%{
u1n(i+1,:) = [ xb_dotf - theta_dotf*cos(thetaf)*(H/2 - S*sin(alpha)) + theta_dotf*sin(thetaf)*(W/2 + S*cos(alpha)) - l1_dotf*cos(alpha)*cos(thetaf) + l1_dotf*sin(alpha)*sin(thetaf) + l1f*theta_dotf*cos(alpha)*sin(thetaf) + l1f*theta_dotf*sin(alpha)*cos(thetaf),...
                  zb_dotf - theta_dotf*cos(thetaf)*(W/2 + S*cos(alpha)) - theta_dotf*sin(thetaf)*(H/2 - S*sin(alpha)) - l1_dotf*cos(alpha)*sin(thetaf) - l1_dotf*sin(alpha)*cos(thetaf) - l1f*theta_dotf*cos(alpha)*cos(thetaf) + l1f*theta_dotf*sin(alpha)*sin(thetaf)...
                  ];
u2n(i+1,:) = [ xb_dotf - theta_dotf*cos(thetaf)*(H/2 - S*sin(alpha)) - theta_dotf*sin(thetaf)*(W/2 + S*cos(alpha)) + l2_dotf*cos(alpha)*cos(thetaf) + l2_dotf*sin(alpha)*sin(thetaf) - l2f*theta_dotf*cos(alpha)*sin(thetaf) + l2f*theta_dotf*sin(alpha)*cos(thetaf),...
                   zb_dotf + theta_dotf*cos(thetaf)*(W/2 + S*cos(alpha)) - theta_dotf*sin(thetaf)*(H/2 - S*sin(alpha)) + l2_dotf*cos(alpha)*sin(thetaf) - l2_dotf*sin(alpha)*cos(thetaf) + l2f*theta_dotf*cos(alpha)*cos(thetaf) + l2f*theta_dotf*sin(alpha)*sin(thetaf)...
                   ]; 
%}
%-----Does Lander contact to the ground?---------           
%if (u1n(i+1,2)*u1n(i,2)<0)
%if (q1n(i,2) <= 0 || (q1n(i,2) < abs(un(i+1,2))*H_Launch_delay*dt + 0.5*g*(H_Launch_delay*dt)^2 && qn(i,5) > 52*pi/180)...
%        || (qn(i, 5)>0 && un(i, 5) > 0))
%if (q1n(i,2) <= 0 || (q1n(i,2) < abs(un(i+1,2))*H_Launch_delay*dt + 0.5*g*(H_Launch_delay*dt)^2 && qn(i,5) > 52*pi/180))
if (q1n(i,2) <= 0)
   if (xmr(i,2)<L_stroke)
      Launchr = true;
   else
      Launchr = false;
   end
end

%if (u2n(i+1,2)*u2n(i,2)<0)
%if(q2n(i,2) <= 0  || (q2n(i,2) < abs(un(i+1,2))*H_Launch_delay*dt + 0.5*g*(H_Launch_delay*dt)^2 && qn(i,5) < -52*pi/180)...
%       || (qn(i, 5)<0 && un(i, 5) < 0))
%if(q2n(i,2) <= 0  || (q2n(i,2) < abs(un(i+1,2))*H_Launch_delay*dt + 0.5*g*(H_Launch_delay*dt)^2 && qn(i,5) < -52*pi/180))
if(q2n(i,2) <= 0)
   if (xml(i,2)<L_stroke)
      Launchl = true;
   else
      Launchl = false;
   end
end
%{
if(Vr(i+1) > 1000)
    Launchr = false;
end
if(Vl(i+1) > 1000)
    Launchl = false;
end
%}
F1_impactn(i) = -k*(qn(i,3)-l1ss)-c*un(i,3);
F2_impactn(i) = -k*(qn(i,4)-l2ss)-c*un(i,4);

end
xmr(end,2) = xmr(end-1,2);
xml(end,2) = xml(end-1,2);
F1max=max(F1_impactn);
F2max=max(F2_impactn);

toc;

if abs(un(end,2))<=1e-2
    flag='stable';
else
    flag='unstable';
end





%% Plot non-linear result
figure(1);
subplot(2,1,1);
plot(t,qn(:,5)); %theta(t)
xlabel('Time [s]','Interpreter','Latex');
ylabel('Theta [rad]','Interpreter','Latex');
title('Theta to Time relation','Interpreter','Latex');
subplot(2,1,2);plot(t,qn(:,2)); %theta(t)
xlabel('Time [s]','Interpreter','Latex');
ylabel('Zb [m]','Interpreter','Latex');
title('Body-height to Time relation','Interpreter','Latex');

figure(2);
subplot(4,1,1);
plot(t,qn(:,2)); %zb(t)
hold on;
plot(t,q1n(:,2)); %xb(t)
plot(t,q2n(:,2)); %l1(t)
plot(t,qn(:,4)); %l2(t)
plot(t,qn(:,5)); %theta(t)
hold off;
L0=legend('$z_{b}$(t) [m]','$x_{b}$(t) [m]','$l_{1}$(t) [m]','$l_{2}$(t) [m]','$Theta$(t) [m]','Interpreter','Latex');
%L0=legend('$Theta$(t) [m]','Interpreter','Latex');
set(L0,'FontSize',10);
grid on;
xlabel('Time [s]','Interpreter','Latex');
ylabel('Response','Interpreter','Latex');
title({'2D Lander Free Fall Landing';['xb(t0)=',num2str(qn(1,1)),'[m],','zb(t0)=',num2str(qn(1,2)),'[m],',...
    '$l_{1}$(t0)=',num2str(qn(1,3)),'[m],','$l_2$(t0)=',num2str(qn(1,4)),'[m],','theta(t0)=',num2str(qn(1,5)),'[m].'];...
    ['  $\dot{xb}$(t0)=',num2str(un(1,1)),'[m/s],','$\dot{zb}$(t0)=',num2str(un(1,2)),'[m/s],',...
    '$\dot{l_{1}}$(t0)=',num2str(un(1,3)),'[m/s].','$\dot{l_{2}}$(t0)=',num2str(un(1,4)),'[m/s],','$\dot{theta}$(t0)=',num2str(un(1,5)),'[m/s].']},'Interpreter','Latex');
%set(findall(gcf,'type','line'),'linewidth',2);
%kfig(2)=figure;

subplot(4,1,2);
plot(t,F1_impactn);hold on;
L0=legend('F1 [N](Impact Force on Footpad1)','Interpreter','Latex');
set(L0,'FontSize',10);
grid on;
xlabel('Time [s]','Interpreter','Latex');
ylabel('Impact Force','Interpreter','Latex');
title(['F1max:',num2str(F1max),'[N]'],'Interpreter','Latex');

subplot(4,1,3);
plot(t,F2_impactn);hold on;
L0=legend('F2 [N](Impact Force on Footpad2)','Interpreter','Latex');
set(L0,'FontSize',10);
grid on;
xlabel('Time [s]','Interpreter','Latex');
ylabel('Impact Force','Interpreter','Latex');
title(['F2max:',num2str(F2max),'[N]'],'Interpreter','Latex');

subplot(4,1,4);
plot(t,un(:,2)); %zb_dot(t)
L0=legend('$\dot{zb}$(t) [m/s]','Interpreter','Latex','Location','Southeast');
set(L0,'FontSize',10);
grid on;
xlabel('Time [s]','Interpreter','Latex');
ylabel('Velocity Response','Interpreter','Latex');
sgtitle(['\bf Lander Landing status : ',flag],'FontSize',10,'Interpreter','Latex');
set(findall(gcf,'type','line'),'linewidth',2);

figure(3);
subplot(2,1,1);
stairs(t,Vr);hold on;
L0=legend('$Vr$(t) [V]','Interpreter','Latex');
set(L0,'FontSize',10);
grid on;
xlabel('Time [s]','Interpreter','Latex');
ylabel('Voltage','Interpreter','Latex');
title('Control Effort','Interpreter','Latex');
set(findall(gcf,'type','stair'),'linewidth',2);

subplot(2,1,2);
plot(t,xmr(:,2));
L0=legend('$xpr$(t) [m]','Interpreter','Latex');
set(L0,'FontSize',10);
grid on;
xlabel('Time [s]','Interpreter','Latex');
ylabel('Right VCM Plate Position [m]','Interpreter','Latex');
title('Right AMEID-VCM Plate Position','Interpreter','Latex');
set(findall(gcf,'type','line'),'linewidth',2);

figure(4);
subplot(2,1,1);
stairs(t,Vl);hold on;
L0=legend('$Vl$(t) [V]','Interpreter','Latex');
set(L0,'FontSize',10);
grid on;
xlabel('Time [s]','Interpreter','Latex');
ylabel('Voltage','Interpreter','Latex');
title('Control Effort','Interpreter','Latex');
set(findall(gcf,'type','stair'),'linewidth',2);

subplot(2,1,2);
plot(t,xml(:,2));
L0=legend('$xpl$(t) [m]','Interpreter','Latex');
set(L0,'FontSize',10);
grid on;
xlabel('Time [s]','Interpreter','Latex');
ylabel('Left VCM Plate Position [m]','Interpreter','Latex');
title('Left AMEID-VCM Plate Position','Interpreter','Latex');
set(findall(gcf,'type','line'),'linewidth',2);

%% -------------Animation----------------
%--------------RGB color-----------------
rgb =[0, 0.4470, 0.7410;...     %new_blue
      0.8500, 0.3250, 0.0980;...%orange
      0.9290, 0.6940, 0.1250;...%soil yellow
      0.4940, 0.1840, 0.5560;...%purple
      0.4660, 0.6740, 0.1880;...%lightly green
      0.3010, 0.7450, 0.9330;...%lightly blue
      0.6350, 0.0780, 0.1840];  %Brown
%----------------------------------------

%----------------------------------------
figure(5);
axis ([-2 2 -0.1 7]); 
axis manual;

for n = 1:40:N
    
    Rc = [cos(qn(n, 5)), -sin(qn(n,5));...
              sin(qn(n,5)), cos(qn(n,5))];
    qb=[qn(n,1);qn(n,2)];
    
    r1c = [-W/2;H/2];
    r2c = [W/2;H/2]; 
    r3c = [-W/2;-H/2];
    r4c = [W/2;-H/2];
    
    r1 = Rc*r1c;
    r2 = Rc*r2c;
    r3 = Rc*r3c;
    r4 = Rc*r4c;
    
    q1 = qb + r1;
    q2 = qb + r2;
    q3 = qb + r3;
    q4 = qb + r4;
    
    plot(qn(n,1),qn(n,2),'ro');
    hold on;
    plot(q1n(n,1),q1n(n,2),'bs');
    plot(q2n(n,1),q2n(n,2),'bs');
    %--------------Lander body-------------
    plot([q1(1) q2(1)],[q1(2) q2(2)],'color',rgb(7,:));
    plot([q2(1) q4(1)],[q2(2) q4(2)],'color',rgb(7,:));
    plot([q4(1) q3(1)],[q4(2) q3(2)],'color',rgb(7,:));
    plot([q3(1) q1(1)],[q3(2) q1(2)],'color',rgb(7,:));
    %--------------------------------------
    plot([q1n(n,1) qn(n,1)],[q1n(n,2) qn(n,2)],'k:');
    plot([q2n(n,1) qn(n,1)],[q2n(n,2) qn(n,2)],'k:');
    %--------------------------------------
    floorx = [-10 10];
    floorz = [0 0];
    line(floorx,floorz,'Color','green','LineStyle','--');
    
    xlabel('x [m]','Interpreter','Latex'); 
    ylabel('z [m]','Interpreter','Latex'); 
     
    
    axis ([-10 10 -0.1 10]); 
    %pbaspect([40 71 1]);
    pbaspect([2 1 1]);
    grid on;
    grid minor;
    title('2D animation','Interpreter','Latex');
    set(findall(gcf,'type','line'),'linewidth',2);
    drawnow 
    
    if n==1
        pause(2);
    end
%       % Capture the plot as an image 
%       frame = getframe(h); 
%       im = frame2im(frame); 
%       [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
%       if n == 1 
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%       else 
%           imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%       end 
hold off;
end

