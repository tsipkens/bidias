function [ll, mm, factor, factorInit, numInitTwomeySteps, maxTwomeySteps, smoothFactor] = InversionParameters

ll = 600; % num points at which to describe the distribution in the first dimension (mobility diameter)
mm = 600; % num points at which to describe the distribution in the second dimension (mass)
factor = 5; % factor to be used in Twomey routine except in special cases
factorInit = 5; % F, initial factor to be used in Twomey routine 
numInitTwomeySteps = 2; % F, number of times to employ initial factor
smoothFactor = 1/50;  %Sf 
maxTwomeySteps = 1000;  %A maximum number of Twomey steps