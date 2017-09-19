# Author: Steven Ruggiero
# Email: Stevemruggiero@gmail.com
# 
# This file contains the function that returns the differentials for
# the deterministic model.
#
# 9-8-16: Changed the T cell equations to match the paper exactly.
#       Changed the other T cell equations to match term thre of dT_1dt
#       Turns out I used the terms for dF_adt for dI_ydt
#
# 9-18-16: This file was branched off from DetModel. This is an attempt to
# resolve differences between the Sud and Wigginton papers' differences
#
# 9-19-16: this file was branched off from WSDetModel. Attempt to recreate 
# The Wigginton model.

def WDetModel(y,t,p):

    M_R, M_I, M_A, T_0, T_1, T_2, I_y, I_12, I_10, I_4, B_I, B_E = y
    
    T_T = T_0 + T_1 + T_2
    
    B_T = B_I + B_E
    
    dM_Rdt = (p['s_m'] +
            p['a_4'] * (M_A + p['w'] * M_I) +
            p['a_21']*M_R*B_T/(B_T+p['c_28']) +
            p['k_4']*M_A * (I_10 / (I_10 + p['s_8'])) +
            -p['k_2']*M_R * (B_E / (B_E + p['c_9'])) +
            -p['k_3']*M_R * (I_y / (I_y + p['s_3'])) * (B_T / (B_T + p['c_8'])) +
             p['u_da']*M_A * p['s_3']/(I_y + p['s_3']) * p['c_8']/(B_T + p['c_8']) +
            -p['u_MR']*M_R)
            
    dM_Idt = (p['k_2']*M_R * (B_E / (B_E + p['c_9']))
             -p['k_17']*M_I * (B_I**p['m'] / (B_I**p['m'] + (p['N']*M_I)**p['m'] + p['e'])) + 
            -p['k_14']*M_I * ((T_T/M_I) / (T_T/M_I + p['c_4'])) *
            (1 - p['p'] * B_I/(B_I + p['N']*M_I + p['e'])) +
            -p['u_MI']*M_I)
        
    dM_Adt = ((p['k_3']*M_R * I_y/(I_y + p['s_3']) * B_T /(B_T + p['c_8'])) +
            -p['k_4']*M_A * (I_10 / (I_10 + p['s_8'])) +
            -p['u_da']*M_A * p['s_3']/(I_y + p['s_3']) * p['c_8']/(B_T + p['c_8']) +
            -p['u_MA']*M_A)
            
    dT_0dt = (p['a_1']*(M_A + p['w']*M_I) +
            p['a_2']*T_0 * M_A/(M_A + p['c_15']) +
            -p['u_T0']*T_0 +
            -p['k_6']*I_12*T_0 * I_y/(I_y + p['f_1']*I_4 + p['s_1']) +
            -p['k_7']*T_0 * I_4 / (I_4 + p['f_2']*I_y + p['s_2']) )
            
    dT_1dt = (p['a_3'] * (M_A + p['w']*M_I) +
            p['k_6']*I_12*T_0 * I_y/(I_y + (p['f_1']*I_4) + p['s_1']) +
            -p['u_T1']*T_1)
            
    dT_2dt = (p['a_3'] * (M_A + p['w']*M_I) +
            p['k_7']*T_0 * I_4 / (I_4 + p['f_2']*I_y + p['s_2']) +
            -p['u_T2']*T_2)    
            
    dI_ydt = (p['s_g'] * B_T/(B_T + p['c_10']) * I_12/(I_12 + p['s_7']) +
            p['a_5']*T_1 * M_A/(M_A + p['c_14']) + 
            p['a_7']*T_0 * I_12/(I_12 + p['f_4']*I_10 + p['s_4']) * M_A/(M_A + p['c_14'])+
            -p['u_iy']*I_y)
            
    dI_12dt = (p['a_8']*M_A +
            p['a_22']*M_I +
            p['a_10']*M_R * I_y/(I_y + p['f_5']*I_10 + p['s_5']) * B_T/(B_T + p['c_18'])
            -p['u_i12']*I_12)

    dI_10dt = ((p['a_13']*M_R * B_T/(B_T + p['c_12']) + p['a_14']*M_A) *
            p['s_6']/(I_10 + p['f_6']*I_y + p['s_6']) +
            p['a_16']*T_1 +
            p['a_17']*T_2 +
            p['a_18']*T_0 * I_12/(I_12 + p['s_9']) +
            -p['u_i10']*I_10)
            
    dI_4dt = (p['a_11']*T_0 +
            p['a_12']*T_2 +
            -p['u_i4']*I_4)
            
    dB_Edt = (p['a_20']*B_E +
            -p['k_15']*M_A*B_E +
            -p['k_18']*M_R*B_E +
            p['u_MI']*B_I +            
            p['k_17']*p['N']*M_I * B_I**p['m'] / (B_I**p['m'] + (p['N']*M_I)**p['m'] + p['e']) +
            -p['k_2'] * p['N']/2 * M_R * B_E/(B_E + p['c_9']))

    dB_Idt = (p['a_19']*B_I * (1 - B_I**p['m']/(B_I**p['m'] + (p['N']*M_I)**p['m'] + p['e'])) +
            -p['k_17']*p['N']*M_I * B_I**p['m']/(B_I**p['m'] + (p['N']*M_I)**p['m'] + p['e']) +
            p['k_2'] * p['N']/2 * M_R * B_E/(B_E + p['c_9']) +
            -p['k_14']*p['N'] * M_I * (T_T / M_I) /((T_T / M_I) + p['c_4']) *
            (1 - p['p']*B_I/(B_I + (p['N']*M_I) + p['e'])) +
            -p['u_MI']*B_I)
            
   
    return [dM_Rdt, dM_Idt, dM_Adt, dT_0dt, dT_1dt, dT_2dt, dI_ydt, dI_12dt, dI_10dt, dI_4dt, dB_Idt, dB_Edt]

           