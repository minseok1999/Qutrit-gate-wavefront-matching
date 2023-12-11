%% Subscript for updating phasescreens 
% using current Beam and BeamBack values.
% Not a function since Beam and BeamBack are huge, and we're not sure
% if matlab uses pointers in functions.

DeltaHologram = zeros(nx,ny);

% Calculate overlap o_kii (BeamBack was conjugated in the
% initialization)
Overlap = BeamBack(:,:,PhScrInd).*Beam(:,:,PhScrInd).*exp(1i*Hologram(:,:,PhScrInd));
% Average phase
AvePhase = mean(angle(Overlap));
% Values for updating the Hologram (For a single mode pair):
DeltaHologram = Overlap.*exp(-1i*AvePhase);

% Calculate the total change required for this hologram:
DeltaHoloSum = -angle(sum(DeltaHologram,3));
%Offset
offset=0.0005;
% New Hologram
Hologram(:,:,PhScrInd) = Hologram(:,:,PhScrInd)+DeltaHoloSum+offset;

% Include phase resolution
Hologram(:,:,PhScrInd) = mod(Hologram(:,:,PhScrInd),2*pi);
% Normalize, scale to max phase value and discretize
Hologram2(:,:,PhScrInd) = floor(Hologram(:,:,PhScrInd)/max(max(Hologram(:,:,PhScrInd)))*MaxPhaseValue);
% Back to radians
Hologram(:,:,PhScrInd) = Hologram2(:,:,PhScrInd)/max(max(Hologram2(:,:,PhScrInd)))*2*pi;

