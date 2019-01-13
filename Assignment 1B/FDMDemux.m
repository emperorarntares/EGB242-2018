function [xdm] = FDMDemux(muxSignal,t,Mag,fshift,Phishift)
%FDMDemux handles an input signal with any number of streams. It removes
%frequency shifts for any number of signals.

for ii = 1:length(fshift)
    xdm(ii,:) = muxSignal*Mag(ii).*cos(2*pi*fshift(ii)*t + Phishift(ii));
end

end

