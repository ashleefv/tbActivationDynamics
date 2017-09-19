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
# 9-22-16: In noticed I droped a - sign on a term in dB_Idt
#
class DetModel:
    
    
    def DifEqs(self,y,t,p):
    
        M_R, M_I, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E = y
        
        B_T = B_I + B_E
        
        epsilon = 0.0
        
        dM_Rdt = (p['sr_m']
                +p['a_4a'] * (M_A + p['w_2'] * M_I) 
                +p['sr_4b'] * (F_a / (F_a + p['f_8']*I_10 + p['s_4b'])) 
                -p['k_2']*M_R * (B_E / (B_E + p['c_9']))
                -p['k_3']*M_R * (I_y / (I_y + p['f_1']*I_4 + p['s_1']))  
                *(B_T + p['B']*F_a) / (B_T + p['B']*F_a + p['c_8']) 
    #            +p['k_4']*M_A * (I_10 / (I_10 + p['s_8'])) 
                -p['u_MR']*M_R)
                
        dM_Idt = (p['k_2']*M_R * (B_E / (B_E + p['c_9']))
                -p['k_17']*M_I * (B_I**2 / (B_I**2 + (p['N']*M_I)**2+epsilon))  
                -p['k_14a']*M_I * ((T_c + p['w_3']*T_1)/M_I) / ((T_c + p['w_3']*T_1)/M_I + p['c_4']) 
                -p['k_14b']*M_I * F_a/(F_a + p['f_9']*I_10 + p['s_4b']) 
                -p['k_52']*M_I * (( (T_c*( T_1/(T_1+p['c_T1']) ) + p['w_1']*T_1 ) /M_I) 
                /(( T_c * (T_1/(T_1+p['c_T1']))+p['w_1']*T_1)/M_I + p['c_52'])) 
                -p['u_MI']*M_I)
                
        dM_Adt= (p['k_3']*M_R * I_y/(I_y + p['f_1']*I_4 + p['s_1'])
                *(B_T + p['B']*F_a)/(B_T + p['B']*F_a + p['c_8']) 
                -p['k_4']*M_A * (I_10 / (I_10 + p['s_8'])) 
                -p['u_MA']*M_A)
        
    #    if t < 150.0:
    #        dM_Adt = 0.0
    #            
    #    if M_A > 10.0:
    #        if dM_Adt > 0.0:
    #            dM_Adt = 0.0
                
        dT_0dt = (p['a_1a']*(M_A + p['w_2']*M_I) 
                +p['sr_1b'] * F_a/(F_a + p['f_8'] * I_10 + p['s_4b2'])
                +p['a_2']*T_0 * M_A/(M_A + p['c_15'])
                #In the paper, term 4 is inconsistant with it's corresponding term (term 3 of dT_1/dt) in the paper.
                -p['k_6']*I_12*T_0 * I_y/(I_y + (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
                #Later equations in the paper suggest the form of the term here is correct.
                #-p['k_6']*I_12*T_0 * I_y/(I_y * (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
                -p['k_7']*T_0 * I_4 / (I_4 + p['f_2']*I_y + p['s_2'])
                -p['u_T0']*T_0)
                
        dT_1dt = (p['a_3a'] * (M_A + p['w_2']*M_I)
                +p['sr_3b'] * F_a/(F_a + p['f_8'] * I_10 + p['s_4b1'])
                #Below is as written in the paper (potential typo)
                +p['k_6']*I_12*T_0 * I_y/(I_y + (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
                #Below agrees with the corresponding term in other equations
                #+p['k_6']*I_12*T_0 * I_y/(I_y*(p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
                -p['u_Ty'] * I_y/(I_y + p['c']) * T_1*M_A 
                -p['u_T1']*T_1)
                
        dT_2dt = (p['a_3a2'] * (M_A + p['w_2']*M_I)
                +p['sr_3b2'] * F_a/(F_a + p['f_8'] * I_10 + p['s_4b1'])
                +p['k_7']*T_0 * I_4 / (I_4 + p['f_2']*I_y + p['s_2'])
                -p['u_T2']*T_2)
    
        dT_80dt = (p['a_1a']*(M_A + p['w_2']*M_I)#eqns say wM_I but no w exists, just w1, w2, w3
                +p['sr_1b'] * F_a/(F_a + p['f_8'] * I_10 + p['s_4b2'])
                +p['a_2']*T_80 * M_A/(M_A + p['c_15'])
                #-p['k_6']*I_12*T_80 * I_y/(I_y * (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1']) +
                -p['k_6']*I_12*T_80 * I_y/(I_y + (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
                -p['u_T80']*T_80)
    
        dT_8dt = (p['m'] * p['a_3ac']*(M_A + p['w_2']*M_I)
                +p['m'] * p['sr_3bc'] * F_a/(F_a + p['f_8'] * I_10 + p['s_4b1'])
                #+p['m'] *p['k_6']*I_12*T_80 * I_y/(I_y * (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
                +p['m'] *p['k_6']*I_12*T_80 * I_y/(I_y + (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
                -p['u_Tcy'] * I_y/(I_y + p['c_c']) * T_8 * M_A
                -p['u_T8']*T_8)
                
        dT_cdt = (p['m'] * p['a_3ac']*(M_A + p['w_2']*M_I)
                +p['m'] * p['sr_3bc'] * F_a/(F_a + p['f_8'] * I_10 + p['s_4b1'])
                #+p['m'] *p['k_6']*I_12*T_80 * I_y/(I_y * (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
                +p['m'] *p['k_6']*I_12*T_80 * I_y/(I_y + (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
                -p['u_Tcy'] * I_y/(I_y + p['c_c']) * T_c * M_A +
                -p['u_Tc']*T_c)#Steve had u_T8 instead of u_Tc
                
        dF_adt = (p['a_30']*M_I 
                +p['a_31']*M_A * (I_y + p['B_2']*B_T)/(I_y + p['B_2']*B_T + (p['f_1']*I_4 + p['f_7']*I_10) + p['s_10'])
                +p['a_32']*T_1
                +p['a_33'] * (T_c + T_8)
                #In the paper, the parameter is written as u_Fa, which isn't included in table III    
                #term5 = -p['u_Fa']*F_a
                #I suspect the subscript is supposed to be u_TNF, as is written below
                -p['u_TNF']*F_a)
                
        dI_ydt = (p['s_g'] * B_T/(B_T + p['c_10']) * I_12/(I_12 + p['s_7']) 
                +p['a_5a']*T_1 * M_A/(M_A + p['c_5a']) 
                +p['a_5b']*T_8 * M_A/(M_A + p['c_5b']) 
                +p['a_5c']*M_I 
                +p['a_7']*T_0 * I_12/(I_12 + p['f_4']*I_10 + p['s_4'])#* M_A/(M_A + 7000.0)
                +p['a_7']*T_80 * I_12/(I_12 + p['f_4']*I_10 + p['s_4'])#* M_A/(M_A + 7000.0) 
                -p['u_iy']*I_y)
        print(dM_Adt)            
        dI_12dt = (p['s_12'] * B_T/(B_T + p['c_230']) 
                +p['a_23']*M_R * B_T/(B_T + p['c_23'])
                +p['a_8']*M_A * p['s']/(p['s'] + I_10)
                -p['u_i12']*I_12)
                
    
    
        dI_10dt = (p['d_7']*M_A * p['s_6']/(I_10 + p['f_6']*I_y + p['s_6']) 
                +p['a_16']*T_1
                +p['a_17']*T_2
                +p['a_18']*(T_8 + T_c)
                -p['u_i10']*I_10)
                
        dI_4dt = (p['a_11']*T_0
                +p['a_12']*T_2 
                -p['u_i4']*I_4)
                
        dB_Idt = (p['a_19']*B_I * (1 - B_I**2/(B_I**2 + (p['N']*M_I)**2+epsilon))
                +p['k_2'] * p['N']/2 * M_R * B_E/(B_E + p['c_9'])
                -p['k_17']*p['N']*M_I * B_I**2/(B_I**2 + (p['N']*M_I)**2+epsilon)
                -p['k_14a']*p['N'] * M_I * ( (T_c + p['w_3']*T_1) / M_I) / 
                ( (T_c + p['w_3']*T_1) / M_I + p['c_4']) 
                -p['k_14b']*p['N']*M_I * F_a / (F_a + p['f_8']*I_10 + p['s_4b'])
                -p['k_52']*p['N']*M_I * (T_c * T_1/(T_1 + p['c_T1']) + p['w_1']*T_1)/M_I /
                ( (T_c * T_1/(T_1 + p['c_T1']) + p['w_1']*T_1 ) /M_I + p['c_52'])
                -p['u_i']*B_I)
    #    if B_I < 1.0:
    #        if dB_Idt < 0.0:
    #            dB_Idt = 0.0
            
                
        dB_Edt = (p['a_20']*B_E +
                p['u_i']*B_I +
                -p['k_15']*M_A*B_E +
                -p['k_18']*M_R*B_E +
                p['k_17']*p['N']*M_I * B_I**2 / (B_I**2 + (p['N']*M_I)**2+epsilon) +
                -p['k_2'] * p['N']/2 * M_R * B_E/(B_E + p['c_9']) +
                (p['k_14a']*p['N'] * p['N_fracc']*M_I * ((T_c + p['w_3']*T_1) / M_I) / 
                ((T_c + p['w_3']*T_1) / M_I + p['c_4'])) +
                p['k_14b']*p['N'] * p['N_fraca']*M_I * F_a / (F_a + p['f_9']*I_10 + p['s_4b']))
        
    #    if B_E < 1.0:
    #        if dB_Edt < 0.0:
    #            dB_Edt = 0.0
    
        return [dM_Rdt, dM_Idt, dM_Adt, dT_0dt, dT_1dt, dT_2dt, dT_80dt, dT_8dt, dT_cdt, dF_adt, dI_ydt, dI_12dt, dI_10dt, dI_4dt, dB_Idt, dB_Edt]

           