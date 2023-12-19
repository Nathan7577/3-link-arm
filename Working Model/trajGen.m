function constVect = trajGen(th_0,th_f,dth_0,dth_f,T)
    mat = [0        0      0    1;
           T^3      T^2    T    1;
           0        0      1    0;
           3*T^2    2*T    1    0];
    
    boundryConds = [th_0; th_f; dth_0; dth_f];
    
    constVect = mat\boundryConds;
end