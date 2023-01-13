function expkdota = expofix(k_point_coord, numerator)
    
    %kdota is represented by k_point_coord / numerator, of size [3, wavevec_totalnum]
    %we need to calculate kdota modulo 1/4, that is kdota = m/4 + frac
    %output is always a column vector
    
    sizeofk = size(k_point_coord);
    k_point_coord = double(k_point_coord(:));
    numerator = double(numerator);
    period_kdota = round(4.0*k_point_coord / numerator);
    fracpart_kdota = (4.0*k_point_coord - period_kdota * numerator) / (4*numerator);
    expkdota = cospi(2.0*fracpart_kdota) + 1i*sinpi(2.0*fracpart_kdota);
    idx = mod(period_kdota,4)==1;
    expkdota(idx) = 1i * expkdota(idx);
    idx = mod(period_kdota,4)==2;
    expkdota(idx) = -expkdota(idx);     
    idx = mod(period_kdota,4)==3;
    expkdota(idx) = -1i * expkdota(idx);
    expkdota = reshape(expkdota, sizeofk);

end