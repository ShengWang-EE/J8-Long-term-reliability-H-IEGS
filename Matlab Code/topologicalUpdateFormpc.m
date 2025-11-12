function [newmpc,pipelineFailue] = topologicalUpdateFormpc(GEresults,mpc,gtd,Info_components,burstPressure_merged, ...
    rupturePressure_merged,delta_merged,wallThickness,dl,M_ng,gamma_ng,R_air,T_gas,psi,kappa)

newmpc = mpc;
pipelineFailue = 0;
[failedPipe,failedSegment,segmentState, failedBranch, failedGsou, failedGen,electricityLoadLevel,gasLoadLevel] = ...
    deal(Info_components(1),Info_components(2),Info_components(3),Info_components(4),...
        Info_components(5),Info_components(6),Info_components(7),Info_components(8));

burstPressure = burstPressure_merged{failedPipe}(segmentState);
rupturePressure = rupturePressure_merged{failedPipe}(segmentState);
delta = delta_merged{failedPipe}(segmentState);

% get the operating pressure, using a simple method
fb = GEresults.Gline(failedPipe,1); tb = GEresults.Gline(failedPipe,2);
pressure_fb = GEresults.Gbus(fb,7); pressure_tb = GEresults.Gbus(tb,7);
pressureForLSE = max(pressure_fb, pressure_tb);

% calculate the limit stat function
LSEsmallLeak = psi*wallThickness(failedPipe)-delta;
LSElargeLeak = kappa * burstPressure - pressureForLSE;
LSErupture = kappa * rupturePressure - pressureForLSE;

% calculate the pressure at the defect location
l_fb = (failedSegment-1)*dl(failedPipe) + 0.5*dl(failedPipe); % the distance between the segment and the "from bus"
length = gtd.Gline(failedPipe,4) / 1000; % km

C_fb_square = mpc.Gline(failedPipe,3)^2 / length * l_fb;
C_tb_square = mpc.Gline(failedPipe,3)^2 - C_fb_square;
pressureForLeak = sqrt( pressure_fb^2 - (pressure_fb^2 - pressure_tb^2)/length*l_fb );



segmentPressure = sqrt( GEresults.Gbus(fb,7) - sign(GEresults.Gline(failedPipe,6)) * ...
    (GEresults.Gline(failedPipe,6))^2 * mpc.Gline(failedPipe,3)^2 );
% update the mpc
if ( (LSEsmallLeak <= 0) && (LSElargeLeak > 0) ) ...% small leak happens
        || ( (LSElargeLeak <= 0) && (LSErupture > 0 ) ) % large leak
    % create a virtual bus
    virtualBusIndex = 21; virtualBusType = 4; 
    fai = 12 + (31-12)/wallThickness(failedPipe) * delta; % diameter of the leak hole
    areaOfLeakHole = pi*(fai/2)^2 / 1e6; % m^3
    virtualGasLoad = areaOfLeakHole * pressureForLeak*1e5 * ( ...
        (M_ng*gamma_ng)/(R_air*T_gas) * ...
        (2/gamma_ng+1)^((gamma_ng+1)/(gamma_ng-1)) )^0.5 /1e6*3600*24; % Mm3/day
    minPressure = 0; maxPressure = min(mpc.Gbus([fb;tb],6));
    virtualBus = [virtualBusIndex, virtualBusType, virtualGasLoad, pressureForLeak,minPressure,maxPressure ];
    newmpc.Gbus = [mpc.Gbus; virtualBus];

    % creating two virtual gas pipelines
    addGline = [fb, 21, sqrt(C_fb_square), 0, mpc.Gline(failedPipe,5);
                21, tb, sqrt(C_tb_square), 0, mpc.Gline(failedPipe,5)];

    newmpc.Gline = [mpc.Gline; addGline];
    % delete the original gas pipe
    newmpc.Gline(failedPipe,3) = 1e-4;
    pipelineFailue = 1;
elseif (LSErupture <= 0) % rupture    (LSElargeLeak <= 0) && 
    newmpc.Gline(failedPipe,3) = 1e-4; % let Cij=0
    pipelineFailue = 2;
elseif (LSEsmallLeak > 0) && (LSElargeLeak > 0) && (LSErupture > 0) % normal
else
    error('other unexpected pipeline failure mode?')
end



end