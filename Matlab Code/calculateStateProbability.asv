function probability = calculateStateProbability(Info_components,rts,Pr_pipe_merged,nGl,nSegmentForPipeline,nPipelineState,t)

[failedPipe,failedSegment,segmentState, failedBranch, failedGsou, failedGen,electricityLoadLevel,gasLoadLevel] = ...
    deal(Info_components(1),Info_components(2),Info_components(3),Info_components(4),...
        Info_components(5),Info_components(6),Info_components(7),Info_components(8));

avaliability_gen = rts.gen(:,1)./(rts.gen(:,1) + rts.gen(:,2));
avaliability_Gsou = rts.Gsou(:,1)./(rts.Gsou(:,1) + rts.Gsou(:,2));
avaliability_branch = rts.branch(:,2)./(rts.branch(:,1) + rts.branch(:,2));

% pipeline segment
probability_pipeline = 1;
for iPipeline = 1:nGl
    for iSegment = 1:nSegmentForPipeline(iPipeline)
        for iSegmentState = 1:nPipelineState(iPipeline)

                if (iPipeline == failedPipe) && (iSegment == failedSegment)% if it is the failed segment
                    probability_pipeline = probability_pipeline * Pr_pipe_merged{iPipeline}(t,segmentState);
                else % 认为其他元件的概率是非故障概率（1-故障概率）还是无所谓状态的概率（1）比较好呢？
                    probability_pipeline = probability_pipeline * sum(Pr_pipe_merged{iPipeline}(t,1));
                    % 认为是低于该故障的segment所在的状态  之前的状态 的概率之和
                end

        end
    end
end

% branch
individualProbability_branch = avaliability_branch;
if failedBranch ~= 0 % there is branch failed
    individualProbability_branch(failedBranch) = 1-avaliability_branch(failedBranch);
end
probability_branch = prod (individualProbability_branch);

% gas source
individualProbability_Gsou = avaliability_Gsou;
if failedGsou ~= 0 % there is branch failed
    individualProbability_Gsou(failedGsou) = 1-avaliability_Gsou(failedGsou);
end
probability_Gsou = prod (individualProbability_Gsou);

% gen
individualProbability_gen = avaliability_gen;
if failedGen ~= 0 % there is branch failed
    individualProbability_gen(failedGen) = 1-avaliability_gen(failedGen);
end
probability_gen = prod (individualProbability_gen);

% electricity load
individualProbability_electricityLoad = rts.electricityLoad(2,electricityLoadLevel);
% gas load
individualProbability_gasLoad = rts.gasLoad(2,gasLoadLevel);

probability = probability_pipeline * probability_branch * probability_Gsou * ...
    probability_gen * individualProbability_electricityLoad * individualProbability_gasLoad;


end