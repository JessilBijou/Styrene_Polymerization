c_A0 = 10;
c_B0 = 1.167;
c_C0 = 0;
V0  = 1;
T0 = 50 + 273.15;
c_Bin = 20;
x0 = zeros(3,1);
x0(1) = c_A0*V0;
x0(2) = x0(1) + c_C0*V0;
x0(3) = V0;
x0(4) = T0;
u0 = zeros(3,1);
u0(2) = 40;
u0(3) = c_Bin;
y0 = fedbatch_OutputFcn(x0,u0);

%% Nonlinear MPC Design to Optimize Batch Operation 
% nonlinear MPC object with 4 states, 3 outputs, 2 manipulated
% variables, and 1 measured disturbance.
nlmpcobj_Plan = nlmpc(4, 3, 'MV', [1,2], 'MD', 3);
Tf = 0.5;
N = 50;
Ts = Tf/N;
nlmpcobj_Plan.Ts = Ts;
nlmpcobj_Plan.PredictionHorizon = N;
nlmpcobj_Plan.ControlHorizon = [7 7 7 7 7 7 7 1];
nlmpcobj_Plan.Model.StateFcn = @(x,u) fedbatch_StateFcnDT(x,u,Ts);
nlmpcobj_Plan.Model.OutputFcn = @(x,u) fedbatch_OutputFcn(x,u);
nlmpcobj_Plan.Model.IsContinuousTime = false;
nlmpcobj_Plan.MV(1).Min = 0;
nlmpcobj_Plan.MV(1).Max = 1;
nlmpcobj_Plan.MV(2).Min = 20;
nlmpcobj_Plan.MV(2).Max = 60;
nlmpcobj_Plan.OV(2).Max = 1.45e5;
nlmpcobj_Plan.OV(3).Max = 1.1;
nlmpcobj_Plan.Optimization.CustomCostFcn = @(X,U,e,data) X(end,1)-X(end,2);
nlmpcobj_Plan.Optimization.ReplaceStandardCost = true;
nlmpcobj_Plan.Optimization.MVInterpolationOrder = 1;
nlmpcobj_Plan.Optimization.SolverOptions.StepTolerance = 1e-8;
validateFcns(nlmpcobj_Plan, x0, u0(1:2), u0(3));

% production of |C| is maximized at the end of the batch process. To do so,
% use the |nlmpcmove| function.
fprintf('\nOptimization started...\n');
[~,~,Info] = nlmpcmove(nlmpcobj_Plan,x0,u0(1:2),zeros(1,3),u0(3));
fprintf('   Expected production (y1) is %g moles.\n',Info.Yopt(end,1));
fprintf('   First order optimality is satisfied (Info.ExitFlag = %i).\n',...
    Info.ExitFlag);
fprintf('Optimization finished...\n');
Nstep = size(Info.Xopt,1) - 1;
t = 0;
X = x0';
t0 = 0;
for i = 1:Nstep
    u_in = [Info.MVopt(i,1:2)'; c_Bin];
    ODEFUN = @(t,x) fedbatch_StateFcn(x, u_in);
    TSPAN = [t0, t0+Ts];
    Y0 = X(end,:)';
    [TOUT,YOUT] = ode15s(ODEFUN,TSPAN,Y0);
    t = [t; TOUT(2:end)];
    X = [X; YOUT(2:end,:)];
    t0 = t0 + Ts;
end
nx = size(X,1);
Y = zeros(nx,3);
for i = 1:nx
    Y(i,:) = fedbatch_OutputFcn(X(i,:)',u_in)';
end
fprintf('\n   Actual Production of C (y1) is %g moles.\n',X(end,2)-X(end,1));
fprintf('   Heat removal rate (y2) satisfies the upper bound.\n');
figure
subplot(2,1,1)
plot(t,Y(:,1),(0:Nstep)*Ts, Info.Yopt(:,1),'*')
axis([0 0.5 0 Y(end,1) + 0.1])
legend({'Actual','Product'},'location','northwest')
xlabel('Time(h)')
ylabel('Product in reactor(mol/m3)')
title('Product formed in reactor (mol/m3)')
subplot(2,1,2)
tTs = (0:Nstep)*Ts;
t(end) = 0.5;
plot(t,Y(:,2),'-',[0 tTs(end)],1.5e5*ones(1,2),'r--')
axis([0 0.5 0.8e5, 1.6e5])
legend({'Heat  removal','Upper Bound'},'location','southwest')
xlabel('Time(h)')
ylabel('Heat removal rate (J.h-1)')
title('Heat removal rate (J.h-1)')
figure
subplot(2,1,1)
stairs(tTs,Info.MVopt(:,1))
xlabel('Time(h)')
ylabel('feed rate (Kg/h)')
title('Feed rate of solvent (Kg/h)')
subplot(2,1,2)
plot(tTs,Info.MVopt(:,2),'*',t,X(:,4)-273.15,'-',...
    [0 0.5],[20 20],'r--',[0 0.5],[50 50],'r--')
axis([0 0.5 20 50])
xlabel('Time(h)')
ylabel('Heat removal rate (J.h-1)')
title('Reactor temperature (Centigrade)')
legend({'Stpt Temp','react Temp'},'location','southeast')
figure
subplot(2,1,1)
c_A = X(:,1)./X(:,3);
c_B = (c_Bin*X(:,3) + X(:,1) + V0*(c_B0 - c_A0 - c_Bin))./X(:,3);
plot(t,[c_A, c_B])
xlabel('Time(h)')
ylabel('concentration(mol/m3)')
title('Concentration of monomer and solvent(mol/m3)')
legend({'concentration of monomer','concentration of solvent'}, 'location', 'west')
subplot(2,1,2)
plot(tTs,Info.Yopt(:,3))
xlabel('Time(h)')
ylabel('Volume(m3)')
title('Liquid volume(m3)')

%% for Tracking the Optimal C Product Trajectory
% To track the optimal trajectory of product |C| calculated above,
% design another nonlinear MPC controller with the same prediction model
% and constraints. However, use the standard quadratic cost and default
% horizons for tracking purposes.
%
% To simplify the control task, assume that the optimal trajectory of the
% |B| feed rate is implemented in the plant and the tracking controller
% considers it to be a measured disturbance. Therefore, the controller uses
% the reactor temperature setpoint as its only manipulated variable to
% track the desired |y1| profile.

%%
% Create the tracking controller.
nlmpcobj_Tracking = nlmpc(4,3,'MV',2,'MD',[1,3]);
nlmpcobj_Tracking.Ts = Ts;
nlmpcobj_Tracking.Model = nlmpcobj_Plan.Model;
nlmpcobj_Tracking.MV = nlmpcobj_Plan.MV(2);
nlmpcobj_Tracking.OV = nlmpcobj_Plan.OV;
nlmpcobj_Tracking.Weights.OutputVariables = [1 0 0];        % track y1 only
nlmpcobj_Tracking.Weights.ManipulatedVariablesRate = 1e-6;  % agressive MV
Cref = Info.Yopt(:,1);
MD = [Info.MVopt(:,1) c_Bin*ones(N+1,1)];
[X1,Y1,MV1,et1] = fedbatch_Track(nlmpcobj_Tracking,x0,u0(2),N,Cref,MD);
fprintf('\nNonlinear MPC: Elapsed time = %g sec. Production of C = %g mol\n',et1,Y1(end,1));
nlmpcobj_Tracking.Optimization.RunAsLinearMPC = 'Adaptive';
[X2,Y2,MV2,et2] = fedbatch_Track(nlmpcobj_Tracking,x0,u0(2),N,Cref,MD);
fprintf('\nAdaptive MPC: Elapsed time = %g sec. Production of C = %g mol\n',et2,Y2(end,1));
nlmpcobj_Tracking.Optimization.RunAsLinearMPC = 'TimeVarying';
[X3,Y3,MV3,et3] = fedbatch_Track(nlmpcobj_Tracking,x0,u0(2),N,Cref,MD);
fprintf('\nTime-varying MPC: Elapsed time = %g sec. Production of C = %g mol\n',et3,Y3(end,1));
figure
plot(Ts*(0:N),[Y1(:,1)])
title('Production of C')
legend({'NLMPC'},'location','northwest')
figure
plot(Ts*(0:N),[Y1(:,2) ])
xlabel('Time(h)')
ylabel('Heat removal rate (J.h-1)')
title('Heat removal rate (J.h-1)')
legend({'NLMPC','Adaptive','TimeVarying','Constraint'},'location','southwest')
