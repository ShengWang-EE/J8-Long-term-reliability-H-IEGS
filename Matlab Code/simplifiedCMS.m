function [simplifiedFlag,electricityLoadCurtailment,gasLoadCurtailment,deviationToACGR] ...
            = simplifiedCMS(Info_components,mpc,smallProbabilityStates,nb,nGb,id,iGd)
simplifiedFlag = 0;
electricityLoadCurtailment = 0; gasLoadCurtailment = 0; deviationToACGR = 0;

[failedPipe,failedSegment,segmentState, failedBranch, failedGsou, failedGen,electricityLoadLevel,gasLoadLevel] = ...
    deal(Info_components(1),Info_components(2),Info_components(3),Info_components(4),...
        Info_components(5),Info_components(6),Info_components(7),Info_components(8));
%%
if smallProbabilityStates == 1% 概率很小
    simplifiedFlag = 1;
    electricityLoadCurtailment = zeros(nb,1);
    gasLoadCurtailment = zeros(nGb,1);
    deviationToACGR = zeros(nGb,1); % !!这里感觉还有争议
elseif failedBranch == 0 && segmentState == 1 % no transmission failure in both systems
    if failedGsou == 0 && failedGen == 0 % no generation or gas source failure
        simplifiedFlag = 1;
        electricityLoadCurtailment = zeros(nb,1);
        gasLoadCurtailment = zeros(nGb,1);
        deviationToACGR = zeros(nGb,1);
    elseif failedGsou == 0 % there is generation failure, but no gas source failure
        simplifiedFlag = 1;

        electricityGeneration = sum(mpc.gen(:,9)); electricityLoad = sum(mpc.bus(:,3));
        gasProduction = sum(mpc.Gsou(:,4)); gasLoad = sum(mpc.Gbus(:,3));
       
        electricityLoadCurtailment_total = (electricityLoad>electricityGeneration).*(electricityLoad-electricityGeneration);
        gasLoadCurtailment_total = (gasLoad>gasProduction) .* (gasLoad-gasProduction);

        electricityLoadCurtailment = zeros(nb,1);
        electricityLoadCurtailment(id) = electricityLoadCurtailment_total / electricityLoad * mpc.bus(id,3);
        gasLoadCurtailment = zeros(nGb,1);
        gasLoadCurtailment(iGd) = gasLoadCurtailment_total / gasLoad * mpc.bus(iGd,3);
        deviationToACGR = zeros(nGb,1);

        % or gas source failure
%         simplifiedFlag = 1;
% 
%         electricityGeneration = sum(mpc.gen(:,9)); electricityLoad = sum(mpc.bus(:,3));
%         gasProduction = sum(mpc.Gsou(:,4)); gasLoad = sum(mpc.Gbus(:,3));
%        
%         electricityLoadCurtailment_total = (electricityLoad>electricityGeneration).*(electricityLoad-electricityGeneration);
%         gasLoadCurtailment_total = (gasLoad>gasProduction) .* (gasLoad-gasProduction);
% 
%         electricityLoadCurtailment = zeros(nb,1);
%         electricityLoadCurtailment(id) = electricityLoadCurtailment_total / electricityLoad * mpc.bus(id,3);
%         gasLoadCurtailment = zeros(nGb,1);
%         gasLoadCurtailment(iGd) = gasLoadCurtailment_total / gasLoad * mpc.bus(iGd,3);
%         deviationToACGR = zeros(nGb,1);
    end
end


end



