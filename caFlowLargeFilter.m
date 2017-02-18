function [uvLarge] = caFlowLargeFilter (uv, threshold)

%discard small magnitude optical flow
magnitudeFlow(:,:,1) = abs(uv(:,:,1)) + abs(uv(:,:,2));
magnitudeFlow(:,:,2) = abs(uv(:,:,1)) + abs(uv(:,:,2));

uvLarge = uv;
uvLarge(magnitudeFlow < threshold) = 0;