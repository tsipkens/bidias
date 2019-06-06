function finalDist = simpleSmooth(distribution,smoothLength)

for i = 1:smoothLength
    distribution(2:end-1) = 0.25*distribution(1:end-2) + 0.25*distribution(3:end) + 0.5*distribution(2:end-1);
    
end

finalDist = distribution;