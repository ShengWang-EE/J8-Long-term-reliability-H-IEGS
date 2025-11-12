function [nodalGasInjectionForSingleComposition,nodalGasInjectionForAllCompositions] = ...
    nodalGasInjection2(PGs,Qptg,omega,Gf,mpc,nGs,nGb,nGl,nPTG,nGasType)
%% nodal injections
% PGs
PGsbus = mpc.Gsou(:,1) ; 
CgsPGs = sparse(PGsbus, (1:nGs)', 1, nGb, nGs); % connection matrix
for r = 1:nGasType
    PGsInbus(:,r) = CgsPGs*(PGs .* mpc.gasCompositionForGasSource(:,r)); % nodal injections of gas source
end
% ptg
PTGbus = mpc.ptg(:,1) ; 
CgsPTG = sparse(PTGbus, (1:nPTG)', 1, nGb, nPTG); % connection matrix
QptgInbusMethane = CgsPTG * Qptg(:,1) ; QptgInbusHydrogen = CgsPTG * Qptg(:,2);
QptgInbus = [QptgInbusMethane, zeros(nGb,3),QptgInbusHydrogen,zeros(nGb,2)];
nodalGasInjectionForSingleComposition = PGsInbus + QptgInbus;
%%
for r = 1:size(Gf,2)
    for i = 1:nGb
        for m = 1:nGl % get all the bus and Gline connected as injection
            if mpc.Gline(m,2) == i % i is the to bus
                % gamma = sign(gasFlow_sum),应该是这个关系
                nodalGasInjectionForSingleComposition(i,r) = nodalGasInjectionForSingleComposition(i,r) + (1+omega(m))/2 * Gf(m,r); 
            elseif mpc.Gline(m,1) == i % i is the from bus
                nodalGasInjectionForSingleComposition(i,r) = nodalGasInjectionForSingleComposition(i,r) + (-1+omega(m))/2 * Gf(m,r);
            end
        end
    end
end
nodalGasInjectionForAllCompositions = sum(nodalGasInjectionForSingleComposition,2);
end