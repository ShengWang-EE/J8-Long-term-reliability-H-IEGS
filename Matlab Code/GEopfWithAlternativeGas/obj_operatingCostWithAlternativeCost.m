function [f,electricityGenerationCost,electricityLoadCurtailmentCost,gasPurchasingCost,gasCurtailmentCost,PTGsubsidy, interchangeabilityCost] = ...
    obj_operatingCostWithAlternativeCost(Pg,PGs,Qptg,LCe,LCg,deviation,mpc,id,iGd)
% unit: $/hour
GCV_ng = 41.04 * 1e6;     % J/m3
CDFe = 1e4; % MW/hour, 大概数值，从jia文章中拿的
% CDFg = CDFe;
CDFg = 1e4 * GCV_ng / 3600 / 24;

penalty_interchangeability = 1e3;

Pg = mpc.baseMVA*Pg;
LCe = mpc.baseMVA*LCe;
%% 
electricityGenerationCost = sum(Pg' * mpc.gencost(:,6));
% electricityLoadCurtailmentCost = sum(LCe) * CDFe;
gasPurchasingCost = sum(PGs' * mpc.Gcost);
% gasCurtailmentCost = sum(LCg) * CDFg;

% 防止一直只切一个节点
gasCurtailmentCost = sum(LCg.*(LCg./mpc.Gbus(iGd,3)/1000+1) * CDFg);
electricityLoadCurtailmentCost = sum(LCe.*(LCe./mpc.bus(id,3)/1000+1) * CDFe);

% subsidy of hydrogen and methane productions
PTGsubsidy = sum(sum(Qptg))  * 1e6/24 * 0.089 * 2.2/ 6.7 * 1;  
% interchangeability cost
interchangeabilityCost = penalty_interchangeability * sum(sum(deviation));
%%

f = electricityGenerationCost + electricityLoadCurtailmentCost + gasPurchasingCost...
    + gasCurtailmentCost - 1*PTGsubsidy + interchangeabilityCost;

end