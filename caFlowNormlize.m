function [uvNorm] = caFlowNormlize(uv)
%normalize optical flow
uvLarge = uv.*100;
varaNorm = sqrt(uvLarge(:,:,1).^2 + uvLarge(:,:,2).^2);
uvNorm(:,:,1) = uvLarge(:,:,1)./varaNorm;
uvNorm(:,:,2) = uvLarge(:,:,2)./varaNorm;
uvNorm(isnan(uvNorm)) = 0;