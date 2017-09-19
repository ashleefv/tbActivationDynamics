# Author: Steven Ruggiero
# Email: Stevemruggiero@gmail.com
# 
# Differential equations for the deterministic model.
# Each differential equation is a seaparte function here for testing and making sure
# the equations are entered correctely, but will be integrated into a single function
# used with the odespy package.
#
#
# 9-8-16: Changed the T cell equations to match the paper exactly.
#
#
#
y =  [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
         #M_R' 
         #M_I'
         #M_A'
         #T_0'
         #T_1
         #T_2
         #T_80
         #T_8
         #T_c
         #F_a
         #I_y
         #I_12
         #I_10
         #I_4
         #B_I
         #B_E
         #B_T #B_T is the sum of B_I and B_E

M_R, M_I, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E = y
    
B_T = B_I + B_E

def dM_Rdt(y,t,p):
#
#
#    
 


    
    term1 = p['sr_m']
    term2 = p['a_4a'] * (M_A + p['w_2'] * M_I)
    term3 = p['sr_4b'] * (F_a / (F_a + p['f_8']*I_10 + p['s_4b']))
    term4 = -p['k_2']*M_R * (B_E / (B_E + p['c_9']))
    term5 = (-p['k_3']*M_R * (I_y / (I_y + p['f_1']*I_4 + p['s_1'])) * 
            ((B_T + p['B']*F_a) / (B_T + p['B']*F_a + p['c_8'])))
    term6 = -p['u_MR']*M_R
    
    return (term1 + term2 + term3 + term4 + term5 + term6, term1, term2, term3, term4, term5, term6)
    

def dM_Idt(y,t,p):
#    
#    
#    
    
    
    
    term1 = p['k_2']*M_R * (B_E / (B_E + p['c_9']))
    term2 = -p['k_17']*M_I * (B_I**2 / (B_I**2 + (p['N']*M_I)**2))
    term3 = -p['k_14a']*M_I * (((T_c + p['w_3']*T_1)/M_I) / ((T_c + p['w_3']*T_1)/M_I + p['c_4']))
    term4 = -p['k_14b']*M_I * F_a/(F_a + p['f_9']*I_10 + p['s_4b'])
    term5 = (-p['k_52']*M_I * (((T_c*(T_1/(T_1+p['c_T1'])) + p['w_1']*T_1)/M_I) /
            ((T_c * (T_1/(T_1+p['c_T1']))+p['w_1']*T_1)/M_I + p['c_52'])))
    term6 = -p['u_MI']*M_I
    
    return (term1 + term2 + term3 + term4 + term5 + term6)
    
def dM_Adt(y,t,p):
#
#
#




    term1 = (p['k_3']*M_R * I_y/(I_y + p['f_1']*I_4 + p['s_1']) *
            (B_T + p['B']*F_a)/(B_T + p['B']*F_a + p['c_8']))
    term2 = -p['k_4']*M_A * (I_10 / (I_10 + p['s_8']))
    term3 = -p['u_MA']*M_A
    
    return (term1 + term2 + term3)
    
    
    
def dT_0dt(y,t,p):
#
#
#




    term1 = p['a_1a']*(M_A + p['w_2']*M_I)
    term2 = p['sr_1b'] * F_a/(F_a + p['f_8'] * I_10 + p['s_4b2'])
    term3 = p['a_2']*T_0 * M_A/(M_A + p['c_15'])
    #In the paper, term 4 is inconsistant with it's corresponding term (term 3 of dT_1/dt).
    #Later equations suggest the form of ther term here is correct.
    term4 = -p['k_6']*I_12*T_0 * I_y/(I_y * (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
    term5 = -p['k_7']*T_0 * I_4 / (I_4 + p['f_2']*I_y + p['s_2'])
    term6 = -p['u_T0']*T_0
    
    return (term1 + term2 + term3 + term4 + term5 + term6, term1, term2, term3, term4, term5, term6)
    
    
    
def dT_1dt(y,t,p):
#
#
#



    term1 = p['a_3a'] * (M_A + p['w_2']*M_I)
    term2 = p['sr_3b'] * F_a/(F_a + p['f_8'] * I_10 + p['s_4b1'])
    #Below is as written in the paper (potential typo)
    term3 = p['k_6']*I_12*T_0 * I_y/(I_y + (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
    #Below agrees with the corresponding term in other equations
    #term3 = p['k_6']*I_12*T_0 * I_y/(I_y*(p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
    term4 = -p['u_Ty'] * I_y/(I_y + p['c']) * T_1*M_A
    term5 = -p['u_T1']*T_1
    
    return (term1 + term2 + term3 + term4 + term5)
    
    
    
def dT_2dt(y,t,p):





    term1 = p['a_3a2'] * (M_A + p['w_2']*M_I)
    term2 = p['sr_3b2'] * F_a/(F_a + p['f_8'] * I_10 + p['s_4b1'])
    term3 = p['k_7']*T_0 * I_4 / (I_4 + p['f_2']*I_y + p['s_2'])
    term4 = -p['u_T2']*T_2
    
    return (term1 + term2 + term3 + term4, term1, term2, term3, term4)
    

def dT_80dt(y,t,p):
    
    
    
    
    #In the paper, the w parameter is missing it's subscript
    #term1 = p['a_1a']*(M_A + p['w']*M_I)
    #Below is what I suspect to be correct
    term1 = p['a_1a']*(M_A + p['w_2']*M_I)
    term2 = p['sr_1b'] * F_a/(F_a + p['f_8'] * I_10 + p['s_4b2'])
    term3 = p['a_2']*T_80 * M_A/(M_A + p['c_15'])
    term4 = -p['k_6']*I_12*T_80 * I_y/(I_y * (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
    term5 = -p['u_T80']*T_80
    
    return (term1 + term2 + term3 + term4 + term5)
    
    
    
def dT_8dt(y,t,p):
    
    
    
    term1 = p['m'] * p['a_3ac']*(M_A + p['w_2']*M_I)
    term2 = p['m'] * p['sr_3bc'] * F_a/(F_a + p['f_8'] * I_10 + p['s_4b1'])
    term3 = p['m'] *p['k_6']*I_12*T_80 * I_y/(I_y * (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1']) 
    term4 = -p['u_Tcy'] * I_y/(I_y + p['c_c']) * T_8 * M_A
    term5 = -p['u_T8']*T_8
    
    return (term1 + term2 + term3 + term4 + term5)
    
def dT_cdt(y,t,p):
    
    term1 = p['m'] * p['a_3ac']*(M_A + p['w_2']*M_I)
    term2 = p['m'] * p['sr_3bc'] * F_a/(F_a + p['f_8'] * I_10 + p['s_4b1'])
    term3 = p['m'] *p['k_6']*I_12*T_80 * I_y/(I_y * (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1']) 
    term4 = -p['u_Tcy'] * I_y/(I_y + p['c_c']) * T_c * M_A
    term5 = -p['u_T8']*T_c
    
    return (term1 + term2 + term3 + term4 + term5)
    
    
def dF_adt(y,t,p):
    
    term1 = p['a_30']*M_I
    term2 = p['a_31']*M_A * (I_y + p['B_2']*B_T)/(I_y + p['B_2']*B_T + (p['f_1']*I_4 + p['f_7']*I_10) + p['s_10'])
    term3 = p['a_32']*T_1
    term4 = p['a_33'] * (T_c + T_8)
    #In the paper, the parameter is written as u_Fa, which isn't included in table III    
    #term5 = -p['u_Fa']*F_a
    #I suspect the subscript is supposed to be u_TNF, as is written below
    term5 = -p['u_TNF']*F_a
    
    return (term1 + term2 + term3 + term4 + term5, term1, term2, term3, term4, term5)
    

def dI_ydt(y,t,p):
    
    term1 = p['s_g'] * B_T/(B_T + p['c_10']) * I_12/(I_12 + p['s_7'])
    term2 = p['a_5a']*T_1 * M_A/(M_A + p['c_5b'])
    term3 = p['a_5b']*T_8 * M_A/(M_A + p['c_5a'])
    term4 = p['a_5c']*M_I
    term5 = p['a_7']*T_0 * I_12/(I_12 + p['f_4']*I_10 + p['s_4'])
    term6 = p['a_7']*T_80 * I_12/(I_12 + p['f_4']*I_10 + p['s_4'])
    term7 = -p['u_iy']*I_y
    
    return (term1 + term2 + term3 + term4 + term5 + term6 + term7)
    
    
    
def dI_12dt(y,t,p):
    
    term1 = p['s_12'] * B_T/(B_T + p['c_230'])
    term2 = p['a_23']*M_R * B_T/(B_T + p['c_23'])
    term3 = p['a_8']*M_A * p['s']/(p['s'] + I_10)
    term4 = -p['u_i12']*I_12
    
    return (term1 + term2 + term3 + term4, term1, term2, term3, term4)
    
def dI_10dt(y,t,p):
    
    term1 = p['d_7']*M_A * p['s_6']/(I_10 + p['f_6']*I_y + p['s_6'])
    term2 = p['a_16']*T_1
    term3 = p['a_17']*T_2
    term4 = p['a_18']*(T_8 + T_c)
    term5 = -p['u_i10']*I_10
    
    return (term1 + term2 + term3 + term4 + term5)
    
def dI_4dt(y,t,p):
    
    term1 = p['a_11']*T_0
    term2 = p['a_12']*T_2
    term3 = -p['u_i4']*I_4
    
    return (term1 + term2 + term3)



# For the next two equations, the paper uses a variable named T_I and paramters C_TI and w_I
# It seems that they've mixed up 1 with I in those cases.
def dB_Idt(y,t,p):

    term1 = p['a_19']*B_I * (1 - B_I**2/(B_I**2 + (p['N']*M_I)**2))
    term2 = p['k_2'] * p['N']/2 * M_R * B_E/(B_E + p['c_9'])
    term3 = -p['k_17']*p['N']*M_I * B_I**2/(B_I**2 + (p['N']*M_I)**2)
    term4 = (-p['k_14a']*p['N'] * M_I * ((T_c + p['w_3']*T_1) / M_I) / 
            ((T_c + p['w_3']*T_1) / M_I + p['c_4']))
    term5 = -p['k_14b']*p['N']*M_I * F_a / (F_a + p['f_8']*I_10 + p['s_4b'])
    #
    term6 = (-p['k_52']*p['N']*M_I * (T_c * T_1/(T_1 + p['c_T1']) + p['w_1']*T_1)/M_I /
            ((T_c * T_1/(T_1 + p['c_T1']) + p['w_1']*T_1)/M_I + p['c_52']))
    term7 = -p['u_i']*B_I
            
    return (term1 + term2 + term3 + term4 + term5 + term6 + term7)
    
    
    
    
    
def dB_Edt(y,t,p):
#
#
#

    
    
    
    term1 = p['a_20']*B_E
    term2 = p['u_i']*B_I
    term3 = -p['k_15']*M_A*B_E
    term4 = -p['k_18']*M_R*B_E
    term5 = p['k_17']*p['N']*M_I * B_I**2 / (B_I**2 + (p['N']*M_I)**2)
    term6 = -p['k_2'] * p['N']/2 * M_R * B_E/(B_E + p['c_9'])
    term7 = (p['k_14a']*p['N'] * p['N_fracc']*M_I * ((T_c + p['w_3']*T_1) / M_I) / 
            ((T_c + p['w_3']*T_1) / M_I + p['c_4']))
    term8 = p['k_14b']*p['N'] * p['N_fraca']*M_I * F_a / (F_a + p['f_9']*I_10 + p['s_4b'])


    return (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8, term1, term2, term3, term4, term5, term6, term7, term8)
    