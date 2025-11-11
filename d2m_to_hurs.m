function hurs = d2m_to_hurs(t2m, d2m) %Temperature in K, Relative Humidity as a % (0-100)
    t2m = t2m - 273.15; %Convert to C
    d2m = d2m - 273.15;
    hurs = 100*(exp((17.625*d2m)./(243.04+d2m))./exp((17.625*t2m)./(243.04+t2m))); % from https://bmcnoldy.earth.miami.edu/Humidity.html
end