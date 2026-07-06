function FHLL = hll_flux(UL,UR,closure,waveMethod)
FL = flux_moment(UL,closure);
FR = flux_moment(UR,closure);

lamL = wave_speeds(UL,closure,waveMethod);
lamR = wave_speeds(UR,closure,waveMethod);

sL = min([lamL(:);lamR(:);0]);
sR = max([lamL(:);lamR(:);0]);

if sL >= 0
    FHLL = FL;
elseif sR <= 0
    FHLL = FR;
else
    FHLL = (sR*FL - sL*FR + sL*sR*(UR-UL))/(sR-sL);
end
end
