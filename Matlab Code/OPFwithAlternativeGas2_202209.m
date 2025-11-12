function [solution, information] = OPFwithAlternativeGas2(mpc,GEresult,contingencyFlag)
%% define named indices into data matrices
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, ...
    TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
% create (read-only) copies of individual fields for convenience\
baseMVA=100;
il = find(mpc.branch(:, RATE_A) ~= 0 & mpc.branch(:, RATE_A) < 1e10);
%% initialization
nGb = size(mpc.Gbus,1);
nGs = size(mpc.Gsou,1);
nGl = size(mpc.Gline,1);
iGd = find(mpc.Gbus(:,3)~=0);
nGd = size(iGd,1);
nb = size(mpc.bus,1);
ng = size(mpc.gen,1);
id = find(mpc.bus(:,3)~=0);
nd = size(id,1);
nPTG = size(mpc.ptg,1);
nGPP = size(find(mpc.gen(:,22)==1),1);
nGasType = size(mpc.gasCompositionForGasSource,2);

omega = sign(GEresult.Gline(:,6)); % direction of the gas flow

refbus = 13;
GCV_ng = 41.04 * 1e6;     % J/m3
M_ng = 17.478 * 1e-3;     % kg/mol
[GCVall, Mall, M_air, fsAll, aAll, R_air, T_stp, Prs_stp, Z_ng, ...
    T_gas, eta_electrolysis,eta_methanation,etaGFU] = initializeParameters();
S_ng = M_ng/M_air;
penalty1 = 100;% penalty factor for SOC relaxations
penalty2 = 1000;% penalty factor for slackvars
penaltyQd = 10000;
penaltyGPP = 10000;
penaltyQ = 10000;
gasComposition_ng = mean(mpc.gasCompositionForGasSource);

gasComposition_original = [
    91.87260333	4.387704675	0.529722883	0.089952937	0.051621302	0.759602627	2.308792207
    89.52850202	5.486821746	0.814856351	0.165418413	0.029984572	0.650791954	3.323624817
    89.52850199	5.486821743	0.814856356	0.165418419	0.029984564	0.65079196	3.323624817
    89.45642031	5.482334802	0.814190034	0.165283193	0.110605057	0.650259806	3.320906883
    91.65999992	3.879999915	0.459999915	0.129999915	0	1.539999915	2.329999915
    91.65999997	3.879999952	0.459999948	0.129999949	0	1.539999947	2.329999953
    92.76569312	2.745670767	0.621558253	0.112365812	0.218893063	0.771823612	2.763994015
    92.18078194	4.319568063	0.429957024	0.02999702	0	0.759924027	2.279772042
    92.18078193	4.319568062	0.429957026	0.029997027	0	0.759924029	2.279772042
    87.62284707	4.105936764	0.408692812	0.028513489	4.944647142	0.722340756	2.167022197
    87.62284703	4.105936737	0.408692789	0.028513469	4.944647099	0.722340733	2.16702217
    87.62284697	4.10593671	0.408692768	0.028513447	4.944647064	0.722340713	2.167022146
    88.56737293	3.780462041	0.376978744	0.027716272	4.481647432	0.759576142	2.006246136
    88.44764581	4.94061048	0.58237266	0.095430216	2.476115568	0.771099468	2.686726085
    88.44764573	4.940610411	0.582372594	0.095430148	2.476115489	0.771099401	2.686726017
    88.44764571	4.940610416	0.582372609	0.095430168	2.47611548	0.77109942	2.686726042
    87.62284711	4.105936732	0.408692762	0.028513451	4.944647087	0.722340697	2.167022151
    87.62284713	4.105936733	0.408692758	0.028513449	4.944647089	0.722340691	2.167022149
    87.62284712	4.105936727	0.408692752	0.028513444	4.944647084	0.722340685	2.167022142
    87.62284712	4.105936737	0.408692766	0.028513459	4.94464707	0.722340697	2.167022152
];
if nGb > 20 % there is virthual bus
    fb = mpc.Gline(20,1);
    tb = mpc.Gline(21,2);
    gasComposition_original_virtualBus = mean([gasComposition_original(fb,:);gasComposition_original(tb,:)]);
    gasComposition_original = [gasComposition_original;gasComposition_original_virtualBus];
end

%% state variables
Prs_square = sdpvar(nGb,1); % bar^2
PGs = sdpvar(nGs,1); % Mm3/day
Qd = sdpvar(nGd,nGasType);% Mm3/day,这里认为是节点消耗的气，加上切负荷后才等于原负荷
Qptg = sdpvar(nPTG,2); % [ methane; hydrogen ] % Mm3/day
Pptg = sdpvar(nPTG,1); % electricity consumption, 1/100 MW
Pg = sdpvar(size(mpc.gen,1),1); % include TPP, GPP and renewable generators, 1/100 MW
Qgpp = sdpvar(nGPP, nGasType); % Mm3/day
Va = sdpvar(nb,1);
gasComposition = sdpvar(nGb,nGasType); 
Gf = sdpvar(nGl,nGasType);% Mm3/day
PHI = sdpvar(nGl,1); %auxiliary variable for gas flow

LCg = sdpvar(nGd,1); % MW
LCe = sdpvar(nd,1); %1/100 MW
gasComposition0 = sdpvar(nGb,nGasType); % the variable in the ACGR
% deviation = sdpvar(nGb,nGasType);
slackVar = sdpvar(nGb,nGasType);
slackQd = sdpvar(nGd,nGasType);
slackGPP = sdpvar(nGPP,nGasType);
slackQ = sdpvar(nGb,nGasType);
%% upper and lower bounds
energyDemand = mpc.Gbus(iGd,3) * GCV_ng; % energy need of these gas bus

Prsmin = mpc.Gbus(:,5); Prsmax = mpc.Gbus(:,6);
PGsmin = mpc.Gsou(:,3); PGsmax = mpc.Gsou(:,4);
Gfmin = -mpc.Gline(:,5); Gfmax = mpc.Gline(:,5); 
LCgmin = 0; LCgmax = energyDemand;
Qptgmin = 0; Qptgmax = mpc.ptg(:,6);
Vamin = -Inf(nb,1); Vamax = Inf(nb,1); 
Vamin(refbus) = 1; Vamax(refbus) = 1;
Pgmin = mpc.gen(:, PMIN) / baseMVA *0; %Pgmin is set to zero
Pgmax = mpc.gen(:, PMAX) / baseMVA;
LCemin = 0; LCemax = mpc.bus(id,PD)/baseMVA;
WImin_normalCondition = 47.2; WImax_normalCondition = 51.41; %MJ/m3
WImin_contingency = 46.5; WImax_contingency = 52.85;
WImin = contingencyFlag*WImin_contingency + (1-contingencyFlag)*WImin_normalCondition;
WImax = contingencyFlag*WImax_contingency + (1-contingencyFlag)*WImax_normalCondition;
ICFmax_normalCondition = 0.48; ICFmin = 0;
ICFmax_contingency = 1.49;
ICFmax = contingencyFlag*ICFmax_contingency + (1-contingencyFlag)*ICFmax_normalCondition;
SImax = 0.6; SImin = 0;
%-------------test------------------------
% WImin = 48; WImax = 50;
% ICFmax = 0.4;

%% constraints
% upper and lower limits
boxCons = [
    Prsmin.^2 <= Prs_square <= Prsmax.^2;
    PGsmin <= PGs <= PGsmax;
        Qd >= 0;    
    LCgmin <= LCg <= LCgmax;
%     Vamin <= Va <= Vamax;
    Pgmin <= Pg <= Pgmax;
    LCemin <= LCe <= LCemax;
    ];

% gas demand cons
QdGasCompositionCons = [];

% for ii = 1:nGd
    QdGasCompositionCons = [
%         -slackQd(ii,:) <= (gasComposition(iGd(ii),:)/100) .* sum(Qd(ii,:),2) - Qd(ii,:) <= slackQd(ii,:); %
        -slackQd <= (gasComposition_original(iGd,:)/100) .* repmat(sum(Qd,2),[1,nGasType]) - Qd <= slackQd;
        ];
% end
gasDemandCons = [
    QdGasCompositionCons;
    slackQd >= 0;
    Qd*GCVall' /1e9 + LCg/1e9*GCV_ng == energyDemand/1e9; % 
%     Qd*1e6/24/3600 * GCVall' /1e6 + LCg*GCV_ng/24/3600 == energyDemand*1e6/24/3600/1e6; % MW
    ]:'gasDemandCons';

% nodal gas flow balance cons
nodalGasFlowBalanceCons = [
    consfcn_nodalGasFlowBalance2(PGs,Qd,Qgpp,Qptg, Gf,mpc,nGasType,nGPP,nGd,iGd) == 0;
    ]:'nodalGasFlowBalanceCons';

% ptg
PTGcons = [
    ( Qptg(:,1) * 1e6/24/3600 *GCVall(1) / eta_methanation + Qptg(:,2) * 1e6/24/3600 * GCVall(5) ...
        ) /1e6 == Pptg*baseMVA / eta_electrolysis; % w
    0 <= Pptg*baseMVA / eta_electrolysis <= Qptgmax /24/3600 * GCVall(5); % 如果全用来制氢，可制备xx
    0 <= Qptg;
    ];

% gpp
Pgpp = Pg(mpc.gfuIndex) * baseMVA;% MW
GPPgasCompositionCons = [];
for ii = 1:nGPP
    GPPgasBus = mpc.GEcon(find(mpc.GEcon(:,2) == mpc.gen(mpc.gfuIndex(ii),1)),1); % find the gas bus for GPP
    GPPgasCompositionCons = [
        GPPgasCompositionCons;
%         (gasComposition(GPPgasBus,:)/100) .* sum(Qgpp(ii,:),2) == Qgpp(ii,:); %
        -slackGPP(ii,:) <= (gasComposition_original(GPPgasBus,:)/100) .* sum(Qgpp(ii,:),2) - Qgpp(ii,:) <= slackGPP(ii,:);
        ];
end
GPPcons = [
    GPPgasCompositionCons;
    Pgpp == Qgpp/24/3600 * GCVall';
    Qgpp >= 0;
    ]:'GPPcons';

% electricity flow
electricityCons = [
    consfcn_electricPowerBalance2(Va,Pg,Pptg,LCe,mpc,id) == 0;
    consfcn_electricBranchFlow2(Va, mpc, il) <= 0;
    Va(refbus) == 1; % 除了slackbus外，其他相角都没约束。一些solver在处理inf的上下限的时候有问题
    ]:'electricityCons';

% security cons 注意这里都是用在ACGR里面的gascomposition来算的
% wobbe index
GCV_nodal0 = gasComposition0/100 * GCVall';
S_nodal0 = gasComposition0/100 * Mall' / M_air;
S_nodal0 = S_ng;
sqrtS0 = 0.5 * (S_nodal0/sqrt(S_ng) + sqrt(S_ng));
WI_nodal0 = GCV_nodal0 ./ sqrtS0 / 1e6; % 如果不能自动转化，那就手动化一下
WobbeIndexCons = [
    WImin * sqrtS0 <= GCV_nodal0/1e6 <= WImax * sqrtS0
    ]:'WobbeIndexCons';
% ICF
gasComposition_np0 = gasComposition0(:,3)/100 + gasComposition0(:,6)/100;
WI_nodal01 = GCV_nodal0 .* ( 1.5*S_ng^(-0.5) - 0.5*S_ng^(-1.5)*S_nodal0 );
ICF_nodal0 = (WI_nodal01 / 1e6 - 50.73+0.03*100*gasComposition_np0)/1.56 - 0.01*100*gasComposition0(:,5)/100;
ICFcons = [
    ICF_nodal0 <= ICFmax;
    ];
% SI
arctanSI_nodal0 = 0.0255*100*gasComposition0(:,3)/100-0.0233*100*gasComposition0(:,6)/100-0.0091*100*gasComposition0(:,5)/100+0.617;
SIcons = [
    arctanSI_nodal0 <= tan(SImax/0.896)
    ];

% SOC reformulation for gas flow
FB = mpc.Gline(:,1); TB = mpc.Gline(:,2);
Gf_sum = sum(Gf,2);
C = mpc.Gline(:,3);
gasFlowSOCcons = [
    (omega-1) .* Gfmax / 2 <= Gf_sum <= (omega+1) .* Gfmax / 2;
    repmat((omega-1) .* Gfmax / 2,[1,nGasType]) <= Gf <= repmat((omega+1) .* Gfmax / 2, [1,nGasType]);
    PHI >= Gf_sum.^2 ./ C.^2 ;
    PHI >= Prs_square(TB) - Prs_square(FB) + (omega + 1) .* (Prsmin(FB).^2 - Prsmax(TB).^2);
    PHI >= Prs_square(FB) - Prs_square(TB) + (omega - 1) .* (Prsmax(FB).^2 - Prsmin(TB).^2);
    PHI <= Prs_square(TB) - Prs_square(FB) + (omega + 1) .* (Prsmax(FB).^2 - Prsmin(TB).^2);
    PHI <= Prs_square(FB) - Prs_square(TB) + (omega - 1) .* (Prsmin(FB).^2 - Prsmax(TB).^2);
%     PHI <= gasFlow_sum0 + 2*gasFlow_sum0 .* (gasFlow_sum - gasFlow_sum0) ./ C.^2 + sigma_PHI;
%     1e4 >= sigma_PHI >= 0;
%     PHI <= 1e4;
    ]:'gasFlowSOCcons';

% gas composition
xi = 0.1;
PGs0 = GEresult.Gsou(:,5);
Qptg0 = [zeros(nPTG,4),GEresult.ptg(:,7),zeros(nPTG,2)];
Gf0 = GEresult.Gline(:,6);
[~,nodalGasInjectionForAllCompositions0] = nodalGasInjection2(PGs0,Qptg0,omega,Gf0,mpc,nGs,nGb,nGl,nPTG,nGasType);
[nodalGasInjectionForSingleComposition,nodalGasInjectionForAllCompositions] = ...
    nodalGasInjection2(PGs,Qptg,omega,Gf,mpc,nGs,nGb,nGl,nPTG,nGasType);

gasCompositionCons = [
%     gasComposition/100 - slackVar <= nodalGasInjectionForSingleComposition./repmat(nodalGasInjectionForAllCompositions0,[1,nGasType]) <= gasComposition/100 + slackVar;
    0 <= slackVar/100 <= 0.1;
    -slackQ <= ...
    repmat(nodalGasInjectionForAllCompositions0,[1,nGasType]).*gasComposition/100 ...
        + repmat(nodalGasInjectionForAllCompositions,[1,nGasType]).* gasComposition_original/100 ...
        - repmat(nodalGasInjectionForAllCompositions0,[1,nGasType]).*gasComposition_original/100 - nodalGasInjectionForSingleComposition ...
        <= slackQ
    slackQ >= 0;
    0 <= gasComposition/100 <= 1;
    0 <= gasComposition0/100 <= 1;
    sum(gasComposition/100,2) == 1;
    sum(gasComposition0/100,2) == 1;
    gasComposition(12,:)/100 == gasComposition(17,:)/100;
    gasComposition0(:,5) <= xi;
%     (1-xi)*GCV_ng/1e6 <= GCV_nodal0/1e6 <= (1+xi)*GCV_ng/1e6;
    ]:'gasCompositionCons';

% summarize all the cons
constraints = [
    boxCons;
    gasDemandCons;
    nodalGasFlowBalanceCons;
    PTGcons;
    GPPcons;
    electricityCons;
    WobbeIndexCons;
    ICFcons;
    SIcons;
    gasFlowSOCcons;
    gasCompositionCons;
    ];
%% solve the problem
nonHydrogen = [1:4,6:7];
penaltyTerm = penalty1 * sum(Gf_sum) + penalty2 * sum(sum(slackVar(:,:)/100)) ...
    + penaltyQd*sum(sum(slackQd(:,:))) + penaltyGPP*sum(sum(slackGPP(:,:))) ...
    + penaltyQ*sum(sum(slackQ));
deviation = (gasComposition-gasComposition0).^2; 
objfcn = obj_operatingCostWithAlternativeCost(Pg,PGs,Qptg,LCe,LCg,deviation,mpc,id,iGd)+ penaltyTerm;
options = sdpsettings('verbose',2,'solver','gurobi', 'debug',1,'usex0',0);
options.gurobi.MIPGap = 5e-2;
information = optimize(constraints, objfcn, options);
%%
%% results
Prs = sqrt(value(Prs_square)); % nodal gas pressure
PGs = value(PGs); % gas production of gas source
Gf = value(Gf); % gas flow in the pipeline
Gf_sum = sum(value(Gf_sum),2);
LCg = value(LCg); % nodal gas load curtailment
Qd = value(Qd); % Mm3/day
Qptg = value(Qptg);
Pptg = value(Pptg); % MW
gasComposition = value(gasComposition); % hydrogen, gas
gasComposition0 = value(gasComposition0);
Va = value(Va); % voltage phase angle
Pg = value(Pg); % electricity generation (including gas fired units)
Pgpp = value(Pgpp);  % MW
Qgpp = value(Qgpp);  % Mm3/day
LCe = value(LCe); % electricity load curtailment
deviation = value(deviation);
deviation_sum = sum(value(deviation),2);
slackVar = value(slackVar);
slackQd = value(slackQd);
slackGPP = value(slackGPP);
slackQ = value(slackQ);
penaltyTerm = value(penaltyTerm);
S_nodal0 = value(S_nodal0);
GCV_nodal0 = value(GCV_nodal0);
WI_nodal0 = value(WI_nodal0);
ICF_nodal0 = value(ICF_nodal0);
SI_nodal0 = atan(value(arctanSI_nodal0))*0.896;
nodalGasInjectionForSingleComposition = value(nodalGasInjectionForSingleComposition);

[objfcn,electricityGenerationCost,electricityLoadCurtailmentCost,gasPurchasingCost,gasCurtailmentCost,PTGsubsidy, interchangeabilityCost] = ...
    obj_operatingCostWithAlternativeCost(Pg,PGs,Qptg,LCe,LCg,deviation,mpc,id,iGd);

[solution.Prs,solution.PGs,solution.Gf,solution.Gf_sum,solution.LCg,solution.Qd,...
    solution.Qptg,solution.Pptg,solution.gasComposition,solution.gasComposition0,...
    solution.Va,solution.Pg,solution.Pgpp,solution.Qgpp,solution.LCe,solution.deviation_sum,...
    solution.objfcn,solution.electricityGenerationCost,solution.electricityLoadCurtailmentCost,...
    solution.gasPurchasingCost,solution.gasCurtailmentCost,solution.PTGsubsidy, solution.interchangeabilityCost,solution.deviation,solution.penaltyTerm] = ...
    deal(Prs,PGs,Gf,Gf_sum,LCg,Qd,Qptg,Pptg,gasComposition,gasComposition0,Va,Pg,Pgpp,Qgpp,LCe,deviation_sum,...
    objfcn,electricityGenerationCost,electricityLoadCurtailmentCost,gasPurchasingCost,gasCurtailmentCost,PTGsubsidy, interchangeabilityCost,deviation,penaltyTerm); 

end