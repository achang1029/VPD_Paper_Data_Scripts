function d2m = hurs_to_d2m(t2m, hurs) % from https://bmcnoldy.earth.miami.edu/Humidity.html
    t2m = t2m - 273.15; %Convert to C

    d2m = 243.04*(log(hurs./100)+((17.625*t2m)./(243.04+t2m)))./(17.625-log(hurs./100)-((17.625*t2m)./(243.04+t2m))); 
    d2m = d2m + 273.15; %Convert back to K
end