function sol = WDetModel(paramFile)
    
	a_1 = [];
	a_2 = [];
	a_3 = [];
	a_4 = [];
	a_5 = [];
	a_7 = [];
	a_8 = [];
	a_10 = [];
	a_11 = [];
	a_12 = [];
	a_13 = [];
	a_14 = [];
	a_16 = [];
	a_17 = [];
	a_18 = [];
	a_19 = [];
	a_20 = [];
	a_21 = [];
	a_22 = [];
	c_4 = [];
	c_8 = [];
	c_9 = [];
	c_10 = [];
	c_12 = [];
	c_14 = [];
	c_15 = [];
	c_18 = [];
	c_28 = [];
	f_1 = [];
	f_2 = [];
	f_4 = [];
	f_5 = [];
	f_6 = [];
	k_2 = [];
	k_3 = [];
	k_4 = [];
	k_6 = [];
	k_7 = [];
	k_14 = [];
	k_15 = [];
	k_17 = [];
	k_18 = [];
	m = [];
	u_i4 = [];
	u_i10 = [];
	u_i12 = [];
	u_iy = [];
	u_da = [];
	u_MA = [];
	u_MI = [];
	u_MR = [];
	u_T0 = [];
	u_T1 = [];
	u_T2 = [];
	N = [];
	p = [];
	s_1 = [];
	s_2 = [];
	s_3 = [];
	s_4 = [];
	s_5 = [];
	s_6 = [];
	s_7 = [];
	s_8 = [];
	s_9 = [];
	s_m = [];
	w = [];
	s_g = [];
	e = [];
    M_R = [];
    M_I = [];
    M_A = [];
    T_0 = [];
    T_1 = [];
    T_2 = [];
    I_y = [];
    I_12 = [];
    I_10 = [];
    I_4 = [];
    B_I = [];
    B_E = [];
    t_0 = [];
    t_f = [];
    t_s = [];
    sourceFile = [];


    load(paramFile);
    
    y = [M_R; M_I; M_A; T_0; T_1; T_2; I_y; I_12; I_10; I_4; B_I; B_E];
    
    sol = ode23t(@WDetDiffs, [t_0, t_f], y);

    function retVector = WDetDiffs(t,y)
        
		
        M_R = y(1);
        M_I = y(2);
        M_A = y(3);
        T_0 = y(4);
        T_1 = y(5);
        T_2 = y(6);
        I_y = y(7);
        I_12 = y(8);
        I_10 = y(9);
        I_4 = y(10);
        B_E = y(11);
        B_I = y(12);

    
        T_T = T_0 + T_1 + T_2;
    
        B_T = B_I + B_E;
    
        dM_Rdt = (s_m + ...
            a_4 * (M_A + w * M_I) + ...
            a_21*M_R*B_T/(B_T+c_28) + ...
            k_4*M_A * (I_10 / (I_10 + s_8)) + ...
            -k_2*M_R * (B_E / (B_E + c_9)) + ...
            -k_3*M_R * (I_y / (I_y + s_3)) * (B_T / (B_T + c_8)) + ...
             u_da*M_A * s_3/(I_y + s_3) * c_8/(B_T + c_8) + ...
            -u_MR*M_R);
            
        dM_Idt = (k_2*M_R * (B_E / (B_E + c_9)) ...
            -k_17*M_I * (B_I^m / (B_I^m + (N*M_I)^m + e)) +  ...
            -k_14*M_I * ((T_T/M_I) / (T_T/M_I + c_4)) * ...
            (1 - p * B_I/(B_I + N*M_I + e)) + ...
            -u_MI*M_I);
        
        dM_Adt = ((k_3*M_R * I_y/(I_y + s_3) * B_T /(B_T + c_8)) + ...
            -k_4*M_A * (I_10 / (I_10 + s_8)) + ...
            -u_da*M_A * s_3/(I_y + s_3) * c_8/(B_T + c_8) + ...
            -u_MA*M_A);
            
        dT_0dt = (a_1*(M_A + w*M_I) + ...
            a_2*T_0 * M_A/(M_A + c_15) + ...
            -u_T0*T_0 + ...
            -k_6*I_12*T_0 * I_y/(I_y + f_1*I_4 + s_1) + ...
            -k_7*T_0 * I_4 / (I_4 + f_2*I_y + s_2));
            
        dT_1dt = (a_3 * (M_A + w*M_I) + ...
            k_6*I_12*T_0 * I_y/(I_y + (f_1*I_4) + s_1) + ...
            -u_T1*T_1);
            
        dT_2dt = (a_3 * (M_A + w*M_I) + ...
            k_7*T_0 * I_4 / (I_4 + f_2*I_y + s_2) + ...
            -u_T2*T_2);   
            
        dI_ydt = (s_g * B_T/(B_T + c_10) * I_12/(I_12 + s_7) + ...
            a_5*T_1 * M_A/(M_A + c_14) + ...
            a_7*T_0 * I_12/(I_12 + f_4*I_10 + s_4) * M_A/(M_A + c_14)+ ...
            -u_iy*I_y);
            
        dI_12dt = (a_8*M_A + ...
            a_22*M_I + ...
            a_10*M_R * I_y/(I_y + f_5*I_10 + s_5) * B_T/(B_T + c_18) ...
            -u_i12*I_12);

        dI_10dt = ((a_13*M_R * B_T/(B_T + c_12) + a_14*M_A) * ...
            s_6/(I_10 + f_6*I_y + s_6) + ...
            a_16*T_1 + ...
            a_17*T_2 + ...
            a_18*T_0 * I_12/(I_12 + s_9) + ...
            -u_i10*I_10);
            
        dI_4dt = (a_11*T_0 + ...
            a_12*T_2 + ...
            -u_i4*I_4);
            
        dB_Edt = (a_20*B_E + ...
            -k_15*M_A*B_E + ...
            -k_18*M_R*B_E + ...
            u_MI*B_I + ...   
            k_17*N*M_I * B_I^m / (B_I^m + (N*M_I)^m + e) + ...
            -k_2 * N/2 * M_R * B_E/(B_E + c_9));

        dB_Idt = (a_19*B_I * (1 - B_I^m/(B_I^m + (N*M_I)^m + e)) + ...
            -k_17*N*M_I * B_I^m/(B_I^m + (N*M_I)^m + e) + ...
            k_2 * N/2 * M_R * B_E/(B_E + c_9) + ...
            -k_14*N * M_I * (T_T / M_I) /((T_T / M_I) + c_4) * ...
            (1 - p*B_I/(B_I + (N*M_I) + e)) + ...
            -u_MI*B_I);
        
        retVector = [dM_Rdt; dM_Idt; dM_Adt; dT_0dt; dT_1dt; dT_2dt; dI_ydt; dI_12dt; dI_10dt; dI_4dt; dB_Edt; dB_Idt;];

    end
end    
    
    