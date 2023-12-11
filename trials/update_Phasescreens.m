%% Subscript for updating phasescreens 
% using current Beam and BeamBack values.
% Not a function since Beam and BeamBack are huge, and we're not sure
% if matlab uses pointers in functions.

DeltaHologram1 = zeros(nx,ny);

% Calculate overlap o_kii (BeamBack was conjugated in the
% initialization)
Overlap1 = BeamBack1(:,:,PhScrInd).*Beam1(:,:,PhScrInd).*exp(1i*Hologram(:,:,PhScrInd));
% Average phase
AvePhase1 = mean(angle(Overlap1));
% Values for updating the Hologram (For a single mode pair):
DeltaHologram1 = Overlap1.*exp(-1i*AvePhase1);


DeltaHologram2 = zeros(nx,ny);

% Calculate overlap o_kii (BeamBack was conjugated in the
% initialization)
Overlap2 = BeamBack2(:,:,PhScrInd).*Beam2(:,:,PhScrInd).*exp(1i*Hologram(:,:,PhScrInd));
% Average phase
AvePhase2 = mean(angle(Overlap2));
% Values for updating the Hologram (For a single mode pair):
DeltaHologram2 = Overlap2.*exp(-1i*AvePhase2);


DeltaHologram3 = zeros(nx,ny);

% Calculate overlap o_kii (BeamBack was conjugated in the
% initialization)
Overlap3 = BeamBack3(:,:,PhScrInd).*Beam3(:,:,PhScrInd).*exp(1i*Hologram(:,:,PhScrInd));
% Average phase
AvePhase3 = mean(angle(Overlap3));
% Values for updating the Hologram (For a single mode pair):
DeltaHologram3 = Overlap3.*exp(-1i*AvePhase3);

%%%%%%%%%%

%%%%%For interpolated optimization

DeltaHologram4 = zeros(nx,ny);

% Calculate overlap o_kii (BeamBack was conjugated in the
% initialization)
Overlap4 = BeamBack4(:,:,PhScrInd).*Beam4(:,:,PhScrInd).*exp(1i*Hologram(:,:,PhScrInd));
% Average phase
AvePhase4 = mean(angle(Overlap4));
% Values for updating the Hologram (For a single mode pair):
DeltaHologram4 = Overlap4.*exp(-1i*AvePhase4);

DeltaHologram5 = zeros(nx,ny);

% Calculate overlap o_kii (BeamBack was conjugated in the
% initialization)
Overlap5 = BeamBack5(:,:,PhScrInd).*Beam5(:,:,PhScrInd).*exp(1i*Hologram(:,:,PhScrInd));
% Average phase
AvePhase5 = mean(angle(Overlap5));
% Values for updating the Hologram (For a single mode pair):
DeltaHologram5 = Overlap5.*exp(-1i*AvePhase5);

DeltaHologram6 = zeros(nx,ny);

% Calculate overlap o_kii (BeamBack was conjugated in the
% initialization)
Overlap6 = BeamBack6(:,:,PhScrInd).*Beam6(:,:,PhScrInd).*exp(1i*Hologram(:,:,PhScrInd));
% Average phase
AvePhase6 = mean(angle(Overlap6));
% Values for updating the Hologram (For a single mode pair):
DeltaHologram6 = Overlap3.*exp(-1i*AvePhase6);


% Calculate the total change required for this hologram:
DeltaHoloSum = -(angle(DeltaHologram1+DeltaHologram2+DeltaHologram3+DeltaHologram4+DeltaHologram5+DeltaHologram6));
% New Hologram
Hologram(:,:,PhScrInd) = Hologram(:,:,PhScrInd)+DeltaHoloSum;

% Include phase resolution
Hologram(:,:,PhScrInd) = mod(Hologram(:,:,PhScrInd),2*pi);
% Normalize, scale to max phase value and discretize
Hologram2(:,:,PhScrInd) = floor(Hologram(:,:,PhScrInd)/max(max(Hologram(:,:,PhScrInd)))*MaxPhaseValue);
% Back to radians
Hologram(:,:,PhScrInd) = Hologram2(:,:,PhScrInd)/max(max(Hologram2(:,:,PhScrInd)))*2*pi;

