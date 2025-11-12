function Info_components = senarioReduction(Info_components,nSystemState)



reducedScenario = [];


for s = 1:nSystemState
    [failedPipe,failedSegment,segmentState, failedBranch, failedGsou, failedGen,electricityLoadLevel,gasLoadLevel] = ...
    deal(Info_components(s,1),Info_components(s,2),Info_components(s,3),Info_components(s,4),...
        Info_components(s,5),Info_components(s,6),Info_components(s,7),Info_components(s,8));
    % 在管道完好，天然气源完好，branch和gen之间只坏一个的情况下，没有切负荷，所以不用算
    if (segmentState==1) && (failedGsou==0) && (failedBranch*failedGen == 0)
        reducedScenario = [reducedScenario;s];
    end
end
        
Info_components(reducedScenario,:) = [];



end