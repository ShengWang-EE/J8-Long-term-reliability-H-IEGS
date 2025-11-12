function [objfcn,electricityGenerationCost,electricityLoadCurtailmentCost,...
    gasPurchasingCost,gasCurtailmentCost] = objfcn_IEGSoperatingCost(Pg,LCe,PGs,Qptg,LCg,mpc,CDFe,CDFg,id,iGd)
Pg = mpc.baseMVA*Pg;
LCe = mpc.baseMVA*LCe;

electricityGenerationCost = sum(Pg' * mpc.gencost(:,6));
% electricityLoadCurtailmentCost = sum(LCe) * CDFe;
gasPurchasingCost = sum(PGs' * mpc.Gcost);
% gasCurtailmentCost = sum(LCg) * CDFg;
PTGsubsidy = sum(sum(Qptg))  *1e6/24 * 0.089 * 2.2/ 6.7;

% 防止一直只切一个节点
gasCurtailmentCost = sum(LCg.*(LCg./mpc.Gbus(iGd,3)/1000+1) * CDFg);
electricityLoadCurtailmentCost = sum(LCe.*(LCe./mpc.bus(id,3)/1000+1) * CDFe);

objfcn = electricityGenerationCost + electricityLoadCurtailmentCost + gasPurchasingCost + gasCurtailmentCost - 1*PTGsubsidy;


end