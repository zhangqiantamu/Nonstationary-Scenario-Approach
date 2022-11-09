function ex_scenario_approach_ppc()
%EX_SCENARIO_APPROACH_PPC  solve the Probabilistic Point Covering via the scenario approach
% 
%   ConvertChanceConstraint (CCC)
%   Copyright (c) 2018-2019
%   by X.Geng
%   Last Edited: July.26.2019
%
%   This file is part of CCC.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%   See https://github.com/xb00dx/ConvertChanceConstraint-ccc for more info.

%% Probabilistic Point Covering (PPC)
% \min_{x,r} r 
%       s.t. P( x-r <= \xi <= x+r ) >= 1-epsilon
%            r >= 0
% \xi is from a standard Gaussian distribution N(0,1)
% this code plots the feasible region and optimal solution

%% Solve PPC via the Scenario Approach
% settings
epsilon = 1e-1; beta = 1e-4;
d = 2; % two decision variable (x,r), also this problem is fully-supported
ops.method = 'scenario approach';
ops.verbose = 0; % debug setting
% prepare data
N = calculate_sample_complexity(d, epsilon, beta, 'Calafiore2006');
% Nonstationary Distribution
shift=0.2;
for i=1:N
wdata_r(:,i)=normrnd(shift*i/N,1+shift*i/N,[1,N]);
xi_scenarios(i)=wdata_r(i,i);
end

% Calculate WS distance
WN = 100000;
testdata=normrnd(shift+1/N,1+shift+1/N,[1,WN]);
for i=1:N
wsd(i) = ws_distance(normrnd(shift*i/N,1+shift*i/N,[1,WN]), testdata, 1);
end
% Make Scenario Robust with Measurement Error Tolerance
rr=2*ones(1,N);
xi_scenarios_up=xi_scenarios+rr;
xi_scenarios_down=xi_scenarios-rr;
xi_robust=[xi_scenarios_up,xi_scenarios_down];

% Traditional Scenario 
sdpvar x r xi
obj = r;
det_constr = [r >= 0];
inner_constr = [x-r <= xi <= x+r];
chance_constr = prob(inner_constr, xi, epsilon, xi_scenarios, ops);
diagnostics = optimize( [chance_constr;det_constr], obj);
disp(diagnostics.info);
assert(diagnostics.problem == 0);
disp( ['optimal x=',num2str(value(x))] );
disp( ['optimal r=',num2str(value(r))] );
% Robust Scenario 
sdpvar x_r r_r xi_r
obj = r_r;
det_constr = [r_r >= 0];
inner_constr = [x_r-r_r <= xi_r <= x_r+r_r];
chance_constr = prob(inner_constr, xi_r, epsilon, xi_robust, ops);
diagnostics = optimize( [chance_constr;det_constr], obj);
disp(diagnostics.info);
assert(diagnostics.problem == 0);
beta_new=solvebeta(N,d,epsilon,wsd,rr);
disp( ['optimal_robust x=',num2str(value(x_r))] );
disp( ['optimal_robust r=',num2str(value(r_r))] );
disp( ['new beta =',num2str(value(beta_new))] );

%% Evaluate the solution
% M = 10000;
% testdata=normrnd(1+1/N,1+1+1/N,[1,M]);
% % evaulate out-of-sample violation probability
% epsilon_test = estimate_violation_probability(inner_constr, xi, testdata, ops);
% disp(['empirical violation probability: ', num2str(epsilon_test)]);
% % find support scenarios
% ops.type = 'convex';
% [sc, sc_indices] = find_support_scenarios(chance_constr, det_constr, obj, xi_scenarios, ops);
% disp([num2str(length(sc_indices)),' support scenarios found: scenario no.', mat2str(sc_indices)]);
% disp(['scenarios: ', mat2str(sc)]);


%% Visualization
% problem description
xi_axis = -4:0.01:5;
ppc_scenario = figure;
Num=5;
for i=1:Num
pdf = normpdf(xi_axis,shift*i/Num,1+shift*i/Num);
plot(xi_axis,pdf,'k-','Color',[1-i/Num,1-i/Num,1]), hold on,
end
% plot(xi_scenarios, zeros(size(xi_scenarios)),'o'), hold on,
% plot(sc, zeros(size(sc)),'rx'), hold on,
% plot( [value(x)-value(r), value(x)+value(r)], [-0.01, -0.01], '^-'), hold on,
% xlabel('\xi'), ylabel('pdf')
% legend('pdf','scenarios','support scenarios','solution')
% title('probabilistic point covering (scenario approach)')
% print(ppc_scenario,'-depsc','-painters','ex_scenario_approach_ppc.eps')
% print(ppc_scenario,'-dsvg','ex_scenario_approach_ppc.svg')

% actual feasible region

% approx feasible region via the scenario approach

end