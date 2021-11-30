%--------------------------------------------------------------------------
clc
clear
%--------------------------------------------------------------------------
set(0, 'DefaultAxesFontName', 'Times');
%% grid model
rng(4)
mpc = loadcase('case33re');

%  Active and reactive loads in the retail market: in [Baren. Wu.]
mpc.bus(:,3) = mpc.bus(:,3)*1; % active power
mpc.bus(:,4) = mpc.bus(:,4)*0; % reactive power

no = length(mpc.bus);   % # of nodes
br = length(mpc.branch);    % # of branch

fn = mpc.branch(:,1);   % from node
tn = mpc.branch(:,2);   % to node

T = graph(fn,tn);   % graph of network topology
Ybus = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);  % make Ybus
Zbus = inv(Ybus(2:end,2:end));    % VSFs
%% network constraints-----change it! true or false

const.activate_Linelimit = true;
const.activate_Voltagelimit = true;
const.activate_Loss = true;

const.Linelimit = [4*ones(1,5) 4*ones(1,6) 1.0*ones(1,21)]';
const.Vmin = 0.95;
const.Vmax = 1.05;

%% Augmented multiplier

rho = 0.02;  % penalty factor

%% For fixed load(nonparticipate p2p energy trade)
%--------------------------------------------------------------------------
Pf = mpc.bus(:,3);  % active power load in the retail market
Qf = mpc.bus(:,4);  % reactive power load in the retail market
pfresult = runpf(mpc,mpoption('verbose',0,'out.all',0));  % compute power flows
Pline_fix = pfresult.branch(:,14);  % active power flows
Qline_fix = pfresult.branch(:,15);  % reactive power flows
V_fix = pfresult.bus(:,8);  % node voltage

%% Agent environment
%--------------------------------------------------------------------------
[agents, sellers, buyers] = market_info33_newfinal(no, rho);   % market information
[centralized_result] = centralized_trading(agents, sellers, buyers, mpc, const);   % centralized_optimization
%--------------------------------------------------------------------------
%% ADMM trading process (Custom setting)
preference = zeros(length(buyers),length(sellers)); % no preferences
preference(4,2) = 0.6; preference(4,3) = 0.3; preference(4,4) = 0;
preference(1,1) = 0.3; preference(1,2) = 0; preference(1,4) = 1.5;
[ADMM.energy, ADMM.sellprice, ADMM.buyprice] = ADMM_trading(agents, sellers, buyers, mpc, const, rho, preference);  

%% ADMM trading process (Unregularized process)
%--------------------------------------------------------------------------
preference = zeros(length(buyers),length(sellers)); % no preferences
const.activate_Linelimit = false;
const.activate_Voltagelimit = false;
const.activate_Loss = false;
[ADMM_unreg.energy, ADMM_unreg.sellprice, ADMM_unreg.buyprice] = ADMM_trading(agents, sellers, buyers, mpc, const, rho, preference);
%--------------------------------------------------------------------------
%% ADMM trading process (Secnario 1 : with preference)
%--------------------------------------------------------------------------
preference = zeros(length(buyers),length(sellers)); % no preferences
preference(4,2) = 0.6; preference(4,3) = 0.3; preference(4,4) = 0;
preference(1,1) = 0.3; preference(1,2) = 0; preference(1,4) = 1.5;
const.activate_Linelimit = true;
const.activate_Voltagelimit = true;
const.activate_Loss = true;
[ADMM_pref.energy, ADMM_pref.sellprice, ADMM_pref.buyprice] = ADMM_trading(agents, sellers, buyers, mpc, const, rho, preference);
%--------------------------------------------------------------------------

ISF=makePTDF(mpc);  % injection shift factor (matpower package)
ISF(:,1) = [];  % remove column index corresponding slack bus

% market layer to physical layer mapping
ADMMtrade_node = agents.As*sum(ADMM.energy,1)'-agents.Ab*sum(ADMM.energy,2);
ADMMtrade_node = ADMMtrade_node/1e3; % unit conversion (kW to MW)

unregADMMtrade_node = agents.As*sum(ADMM_unreg.energy,1)'-agents.Ab*sum(ADMM_unreg.energy,2);
unregADMMtrade_node = unregADMMtrade_node/1e3; % unit conversion (kW to MW)

mpc_ADMMtrade = mpc;
mpc_ADMMtrade.bus(:,3) = mpc_ADMMtrade.bus(:,3) - [0;ADMMtrade_node];
mpc_ADMMtrade = runpf(mpc_ADMMtrade);

mpc_unregADMMtrade = mpc;
mpc_unregADMMtrade.bus(:,3) = mpc_unregADMMtrade.bus(:,3) - [0;unregADMMtrade_node];
mpc_unregADMMtrade = runpf(mpc_unregADMMtrade);
% Assess voltage and line flow

ADMM.VM = mpc_ADMMtrade.bus(:,8);
ADMM.Pline = mpc_ADMMtrade.branch(:,14);

ADMM_unreg.VM = mpc_unregADMMtrade.bus(:,8);
ADMM_unreg.Pline = mpc_unregADMMtrade.branch(:,14);

% Assess total trading energy amount
ADMM.sell = sum(ADMM.energy,1)';
ADMM.buy = sum(ADMM.energy,2);

ADMM_unreg.sell = sum(ADMM_unreg.energy,1)';
ADMM_unreg.buy = sum(ADMM_unreg.energy,2);

main_result_graph(T, no, br, agents, sellers, buyers, ADMM, ADMM_unreg, ADMM_pref, centralized_result, mpc, const, preference);

exactloss = sum(pfresult.branch(:,14)+pfresult.branch(:,16));
mpc.bus(:,3) = mpc.bus(:,3)-[0;ADMMtrade_node];
pfresult = runpf(mpc,mpoption('verbose',0,'out.all',0));
exactloss = sum(pfresult.branch(:,14)+pfresult.branch(:,16));



