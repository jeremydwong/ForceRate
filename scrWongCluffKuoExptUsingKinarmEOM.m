% scrWongCluffKuoExptUsingKinarmEOM
% Goal: generate a coefficient for force-rate that produces minimum jerk
% trajectories and power costs that are the same as the experiment.
% Idea: using Kinarm EOM, compute the joint torques that minimize the net
% energetic cost of reaching, while generating power consumption that is
% very similar to our empirical data.
% this means we want:
% average power to be
%
% if the cForceRateRate is at 1e-4,
% HISTORICAL NOTES ABOUT THIS CODE.
% discovered that the model had damping at the joints in old EOM from DAK.
% discovered that this really helps the system, added them back in and turn
% them down during optimization.
%
% Answer: 3.5e-6 puts force rate cost in the right ballpark. 
clear all;

cForceRateRate = 3.5e-6; %3.5e-6 puts force rate cost in the ballpark. 
cReward = 1;
cRestingRate = 1;
cObj=.0001;


warning off

% load params_TOMLAB

% FROM JDW CLUFF KUO, the target locations and durations. 
fi_targets = [62.0000   87.0000
  66.3376   82.6624
  68.6401   80.3599
  70.0806   78.9194
  70.9822   78.0178].*(pi/180);
sTimeExptHalfCycle =[0.8571    0.6451    0.5173    0.4286    0.3681];

% LOOP THROUGH ALL TARGETS
for iTGT = 1:length(fi_targets)
  
  fi_start  = [fi_targets(iTGT,1),fi_targets(iTGT,1)+5/180*pi];
  fi_target = [fi_targets(iTGT,2),fi_targets(iTGT,2)+5/180*pi];
  %time = durs(TGT);
  
  %% calculating initial states etc...
  F=[0 0]; %no external forces acting on hand
  %[stim_start, q_start, gamma_start, lce_start, Fdes, Ttot] = calc_stable_start2(fi_start,F,grav_angle,l);
  stim_start=[0,0];
  disp('calculated equilibrium start... now running OC problem')
  dampingLevel = [50:-10:0];
  for iDamping = 1:length(dampingLevel) %LOOP THROUGH 3 DAMPING LEVELS, which seems to help. 
    if iDamping==1
      tic
    end
    disp(['Round: ',num2str(iDamping)])
    
    %%%setup the tomlab problem by defining the
    %%%1 TIME
    %%%2 PHASE
    %%%3 CONTROLS
    %%%4 STATES
    toms t %
    tend = sTimeExptHalfCycle(iTGT);
    steps=10;
    p = tomPhase('p', t, 0, tend, steps, [], 'spline3');
    setPhase(p);
    fi = tomState('fi',1,2);        % segment angles
    fid = tomState('fid',1,2);      % segment angular velocities
    stim = tomControl('stim',1,2);  % jer: torque
    %%%/setup the tomlab problem by defining the
    
    %%% initial guess
    if iDamping==1
      %initial guess
      x0 = {
        %tend == 1 % WARNING: if T IS AN OPT VARIABLE, then THIS NEEDS TO BE FIRST IN LIST. propt documentation.
        icollocate({
        fi  == (fi_start-fi_target)/2;%vec(interp1([0 time],[fi_start; fi_target],t))'
        fid  == (fi_target-fi_start)/1;%jdw hack
        })
        collocate({
        %stim  == stim_start %relative muscle activation
        %         stim   == [vec(interp1([0 time],[4.5; 4.51],t/2))',vec(interp1([0 time],[-4.5; -4.51],t/2+1))']
        })
        };
      %%% otherwise warmstart
    else
      x0 = {
        %tend == tendopt %%here we are not solving for this!!
        icollocate({
        fi == fiopt
        fid == fidopt
        })
        
        collocate({
        stim == stimopt
        })
        };
    end
    
    %%% Boundary constraints
    cbnd = {
      initial({
      fi == fi_start;
      fid == 0;
      stim == stim_start;
      })
      final({fi == fi_target;
      fid == 0;
      dot(fid)==0;
      stim == stim_start;
      })
      };
    
    %%% Box constraints
    cbox = {
      %       mcollocate(0 <= fi(1) <= 3*pi/2)
      %       mcollocate(-pi/4 <= fi(2)-fi(1) <= pi)
      %mcollocate(-100 <= fid <= 100)
      
      %0.005 <=  collocate(stim)  <= 1
      %0 <=tend <=10
      };
    %%%
    % ODEs and path constraints via virtual power
    
    load paramsKinarmAlex.mat
    m1 = P.L1_M;
    m2 = P.L2_M;
    m3 = P.L3_M;
    m4 = P.L4_M;
    
    I1 = P.L1_I;
    I2 = P.L2_I;
    I3 = P.L3_I;
    I4 = P.L4_I;
    
    
    l1 = P.L1_L;
    l3 = P.L3_L;
    cx1 = P.L1_C_AXIAL;
    ca1 = P.L1_C_ANTERIOR;
    cx2 = P.L2_C_AXIAL;
    ca2 = P.L2_C_ANTERIOR;
    cx3 = P.L3_C_AXIAL;
    ca3 = P.L3_C_ANTERIOR;
    cx4 = P.L4_C_AXIAL;
    ca4 = P.L4_C_ANTERIOR;
    Q25 = P.L2_L5_ANGLE;
    
    M11 = m1*ca1^2 + m4*ca4^2 + m1*cx1^2 + m4*cx4^2 + m2*l1^2 + I1 + I4;
    M12 = cx2*l1*m2*cos(fi(1) - fi(2)) - ca4*l3*m4*sin(Q25 + fi(1) - fi(2)) + ca2*l1*m2*sin(fi(1) - fi(2)) + cx4*l3*m4*cos(Q25 + fi(1) - fi(2));
    M21 = M12; %Kane/Lagrange derivations produce a symmetric mass matrix.
    M22 = m2*ca2^2 + m3*ca3^2 + m2*cx2^2 + m3*cx3^2 + m4*l3^2 + I2 + I3;
    
    F1 = ca2*l1*m2*fid(2)^2*cos(fi(1) - fi(2)) - cx2*l1*m2*fid(2)^2*sin(fi(1) - fi(2)) - ca4*l3*m4*fid(2)^2*cos(Q25 + fi(1) - fi(2)) - cx4*l3*m4*fid(2)^2*sin(Q25 + fi(1) - fi(2));
    F2 = cx2*l1*m2*fid(1)^2*sin(fi(1) - fi(2)) - ca2*l1*m2*fid(1)^2*cos(fi(1) - fi(2)) + ca4*l3*m4*fid(1)^2*cos(Q25 + fi(1) - fi(2)) + cx4*l3*m4*fid(1)^2*sin(Q25 + fi(1) - fi(2));
    
    MassMat = [M11,M12;M21,M22];
    gradualDamping = [dampingLevel(iDamping)*fid(1)-dampingLevel(iDamping)*fid(2);dampingLevel(iDamping)*fid(2)];
    ceq = collocate({
      dot(fi') == fid'
      MassMat*dot(fid') == [F1;F2] + gradualDamping + [stim(1) - stim(2);stim(2)];
      });
    
    %%% Objective
    
    
    %%% Work cost
    powerSho = stim(1)*fid(1) - stim(2)*fid(1);
    powerElb = stim(2)*fid(2);
    posPowerSho = integrate((powerSho+abs(powerSho))/2);
    posPowerElb = integrate((powerElb+abs(powerElb))/2);
    costPositiveWork = posPowerSho + posPowerElb;
    kMechanicalWorkEfficiencyMargaria = 1/0.25;
    costPositiveWork = costPositiveWork*kMechanicalWorkEfficiencyMargaria;
    
    %%% Force rate cost
    nnshift = 0.00001;
    %       costForceRateLin = cForceRateRate*sum(integrate(sqrt((dot(dot(stim'))).^2+nnshift)));
    costForceRateQuad = cForceRateRate*sum(integrate((dot(dot(stim'))).^2+nnshift));
    costForceRate = costForceRateQuad;
    
    %%% Resting cost
    restingMetRate = 100;%100 Watts.
    costResting = sqrt(tend^2+0.00001)*restingMetRate;
    nameOpt = 'WorkForceRateMin'
    
    
    %%% Reward cost
    reward = 0;%cReward * sqrt(tend^2+0.00001); %avoiding tend<0?
    objective = cObj * (costForceRate + costPositiveWork + reward);
    
    nameSim='Kinarm';
    %% Solve the problem
    options = struct;
    options.name = [nameSim,' ',nameOpt,'.m'];
    %       options.scale ='auto';
    
    %%% optimization setup
    options.PriLev = 2;
    options.Prob.SOL.optPar(35) = 100000; %iterations limit
    options.Prob.SOL.optPar(30) = 200000; %major iterations limit
    Prob = sym2prob(objective, {cbox, cbnd, ceq}, x0, options);
    result = tomRun('snopt',Prob,1);
    solution = getSolution(result);
    fiopt = subs(fi, solution);
    fidopt = subs(fid, solution);
    stimopt = subs(stim, solution);
    
  end %end dampingLevel loop
  
  %%%add to solution vector
  solutions(iTGT) = solution;
  %%%get values
  ECostForceRate(iTGT) = subs(costForceRate,solution);
  ECostPositiveWork(iTGT) = subs(costPositiveWork,solution);
  ECostResting(iTGT) = subs(costResting,solution);
  EReward(iTGT) = subs(reward,solution);
  
  EdotCostForceRate(iTGT) = subs(costForceRate,solution)./tend;
  EdotCostPositiveWork(iTGT) = subs(costPositiveWork,solution)./tend;
  EdotCostResting(iTGT) = subs(costResting,solution)./tend;
  EdotReward(iTGT) = subs(reward,solution)./tend;
  times(iTGT) = tend;
  fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf('TARGET %i\n',iTGT);
  fprintf('Fdotdot power = %.4f \n',EdotCostForceRate(iTGT));
  fprintf('Avg positive power = %.4f\n',EdotCostPositiveWork(iTGT));
  fprintf('Resting power = %.4f\n',EdotCostResting(iTGT));
  fprintf('Reward Power = %.4f\n',EdotReward(iTGT));
  
  s=[fi(1) fid(1) fi(2) fid(2)];
  ang = 0;%grav_angle*180/pi;
  f_k=result.f_k;
  Exit=[result.ExitFlag result.Inform];
  
  %%
  t_plot = linspace(0,subs(tend,solution),451)'; % Time vector where we want our trajectory evaluated.
  %u_plot = collocate(t);
  % figure(1);
  s_plot = atPoints(t_plot,subs(s,solution));
  stim_u = atPoints(t_plot,subs(stim,solution));
  fiplot=s_plot(:,[1 3]);
  fipplot=s_plot(:,[2 4]);
  % [e,h]=showmov3(fi,l,s_plot,0);
  % plotvs(h)
  % axis equal
  figure(8); hold on;
  plot(t_plot,fipplot(:,1),'color',matlabColor(iTGT));
  ylabel('Velocity (m/s)'); xlabel('Time (s)');
  figure(9); hold on;
  plot(t_plot,stim_u(:,1),'color',matlabColor(iTGT));
  plot(t_plot,stim_u(:,2),'color',matlabColor(iTGT));
  ylabel('Force (N)');xlabel('Time (s)');
  tor = atPoints(t_plot,subs(stim,solution));
  t=t_plot;
  state = [fiplot fipplot];
  figure(10); hold on;
  plot(t_plot,fiplot,'color',matlabColor(iTGT));
  xlabel('Position (m)');xlabel('Time (s)');
  
end
%% plot forcerate+work
figure(11);
freqs = 1./sTimeExptHalfCycle;
plot(freqs,EdotCostPositiveWork,'marker','.','markersize',20);hold on;
plot(freqs,EdotCostForceRate,'marker','.','markersize',20);
plot(freqs,EdotCostForceRate+EdotCostPositiveWork,'marker','.','markersize',20);
legend('work','frate','net');

ylabel('Power (W)');
xlabel('Frequency (Hz)');

cd optim_data
eval(['save ',nameOpt,nameSim,' ECostResting ECostPositiveWork ECostForceRate EReward EdotCostResting EdotCostPositiveWork EdotCostForceRate EdotReward times stim t tend steps s fi fid objective f_k Exit solution'])
cd ..

toc
%%
cd optim_data;
figure(8);
savefig([nameOpt,nameSim,'_dxdt']);
figure(9);
savefig([nameOpt,nameSim,'_u']);
figure(10);
savefig([nameOpt,nameSim,'_x']);
figure(11);
savefig([nameOpt,nameSim,'_energybreakdown']);

cd ../
close all;