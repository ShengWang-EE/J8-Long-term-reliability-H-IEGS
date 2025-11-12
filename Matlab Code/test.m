
t=1;s=10;
stateProbability(t,s) = calculateStateProbability(Info_components(s,:),rts,Pr_pipe_merged,nGl,nSegmentForPipeline,nPipelineState,t);
%%
[EDNS,LOLP] = deal(zeros(nb,T)); [EGNS,LOGP,EGID,GIDP] = deal(zeros(nGb,T));
for t = 1:T
    stateCounter = 0;
    for s = 1:nSystemState
        if smallProbabilityStates(s) == 0
            stateCounter = stateCounter + 1;
            EDNS(:,t) = EDNS(:,t) + electricityLoadCurtailment(:,stateCounter) .* stateProbability(t,stateCounter);
            LOLP(:,t) = LOLP(:,t) + ( electricityLoadCurtailment(:,stateCounter) > 1e-4 ) * stateProbability(t,stateCounter);
            EGNS(:,t) = EGNS(:,t) + gasLoadCurtailment(:,stateCounter) * stateProbability(t,stateCounter);
            LOGP(:,t) = LOGP(:,t) + ( gasLoadCurtailment(:,stateCounter) > 1e-4 ) * stateProbability(t,stateCounter);
            EGID(:,t) = EGID(:,t) + deviationToACGR(:,stateCounter) * stateProbability(t,stateCounter);
            GIDP(:,t) = GIDP(:,t) + ( deviationToACGR(:,stateCounter) > 1e-4 ) * stateProbability(t,stateCounter);
        end
    end
end
%% 计算nodal gas security
GCV_nodal0 = gasComposition/100 * GCVall';
S_nodal0 = gasComposition/100 * Mall' / M_air;
S_nodal0 = S_ng;
sqrtS0 = 0.5 * (S_nodal0/sqrt(S_ng) + sqrt(S_ng));
WI_nodal0 = GCV_nodal0 ./ sqrtS0 / 1e6; % 如果不能自动转化，那就手动化一下

gasComposition_np0 = gasComposition(:,3)/100 + gasComposition(:,6)/100;
WI_nodal01 = GCV_nodal0 .* ( 1.5*S_ng^(-0.5) - 0.5*S_ng^(-1.5)*S_nodal0 );
ICF_nodal0 = (WI_nodal01 / 1e6 - 50.73+0.03*100*gasComposition_np0)/1.56 - 0.01*100*gasComposition(:,5)/100;

arctanSI_nodal0 = 0.0255*100*gasComposition(:,3)/100-0.0233*100*gasComposition(:,6)/100-0.0091*100*gasComposition(:,5)/100+0.617;
SI_nodal0 = atan(value(arctanSI_nodal0))*0.896;
%%
computationTime = 0;
for i = 1:size(yalmipTime,1)
    computationTime = computationTime + yalmipTime{i};
end
%%
[EDNS,LOLP] = deal(zeros(nb,T)); [EGNS,LOGP,EGID,GIDP] = deal(zeros(nGb,T));
for t = 1:20
    for s = 1:nReducedScenario
            EDNS(:,t) = EDNS(:,t) + electricityLoadCurtailmentMatrix(:,s) .* stateProbabilityReduced(t,s);
            LOLP(:,t) = LOLP(:,t) + ( electricityLoadCurtailmentMatrix(:,s) > 1e-4 ) * stateProbabilityReduced(t,s)/24;
            EGNS(:,t) = EGNS(:,t) + gasLoadCurtailmentMatrix(:,s) * stateProbabilityReduced(t,s);
            LOGP(:,t) = LOGP(:,t) + ( gasLoadCurtailmentMatrix(:,s) > 1e-4 ) * stateProbabilityReduced(t,s)/24;
            EGID(:,t) = EGID(:,t) + sqrt(deviationToACGRmatrix(:,s)) * stateProbabilityReduced(t,s);
            GIDP(:,t) = GIDP(:,t) + ( deviationToACGRmatrix(:,s) > 1e-4 ) * stateProbabilityReduced(t,s)/24;
    end
end
%% 涉及到LOLP之类概率的都不改了
NJ = size(a,2);
b = a;
for j = 2:NJ
    b(:,j) = ( b(:,j) - b(:,1) ) / 24 + b(:,1);
end
