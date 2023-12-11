function LGout = GenModesLG(ModeTypes, waist, Rad, Angle)
%% This function calculates Laguerre-Gauss modes 
% The created mode indices should be listed in the matrix "ModeTypes", 
% where the first column should include the OAM values and the second 
% column should include the corresponding p-values. 
% - waist should be the wanted beam waist used in calculation
% - Rad should be the wanted 2D radial coordinate grid
% - Angle contains the corresponding azimuthal coordinates
%
% Output is a 3D matrix containing all of the input mode structures stacked
% together in the 3rd dimension of the matrix


ModeNum=size(ModeTypes);
GridDim = size(Rad);
LGout = zeros(GridDim(1), GridDim(2), ModeNum(1));

for ModeInd = 1:ModeNum(1) % Loop through all of the listed modes
         % Extract OAM and p-value
         l = ModeTypes(ModeInd,1); 
         p = ModeTypes(ModeInd,2);
        
         % Calculate the Mode using the Laguerre polynomial function:
         LGTemp = ((sqrt(2)*Rad/waist).^(abs(l))) .* exp(-(Rad.*Rad)/(waist)^2)...
             .* polyval(LaguerreGen(p, abs(l)), 2*(Rad.*Rad)/(waist^2))...
             .* exp(1i*l*Angle);
         
         % Normalize the field:
         NormA = LGTemp.*conj(LGTemp);
         SumOfA = sum(sum(NormA));
         A = sqrt(real(SumOfA));
         
         % Store the value to the output:
         LGout(:, :, ModeInd) = LGTemp/A;          
end
