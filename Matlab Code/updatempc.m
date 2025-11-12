function newmpc1 = updatempc(mpc,Info_components,rts)

[failedPipe,failedSegment,segmentState, failedBranch, failedGsou, failedGen,electricityLoadLevel,gasLoadLevel] = ...
    deal(Info_components(1),Info_components(2),Info_components(3),Info_components(4),...
        Info_components(5),Info_components(6),Info_components(7),Info_components(8));

% pipeline is handled after running the OPF, not here
newmpc1 = mpc;
if failedBranch ~= 0
    % 注意故障的线路不能让capacity=0，否则会造成线路两端电压相等的约束存在
    newmpc1.branch(failedBranch,4) = 1e6; % 应该让电导=0
end
if failedGsou ~= 0
    newmpc1.Gsou(failedGsou,[3,4]) = 0;
end
if failedGen ~= 0
    newmpc1.gen(failedGen,[4:5,9:10]) = 0;
end
newmpc1.bus(:,3:4) = mpc.bus(:,3:4) * rts.electricityLoad(1,electricityLoadLevel);
newmpc1.Gbus(:,3) = mpc.Gbus(:,3) * rts.gasLoad(1,gasLoadLevel);

end