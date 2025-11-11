function y = d2m_to_svp(x) %Temperature in K, Pressure in Pa
    y = 610.8 * exp(17.27 * (x - 273.15) ./ (x - 35.85));
end