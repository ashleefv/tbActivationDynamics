# Author: Steven Ruggiero
# Email: Stevemruggiero@gmail.com
# 
# This file contains the function that returns the differentials for
# the deterministic model, including the colagen model.

import numpy as np
import sys

#ColDetModel is the class that contains the differential equations of the model
#Each differential equation is a separate method which is called by DiffEqs
#self.parameterAdjustFuncs is a dictionary that uses a string representing a parameter
#as the key, and each entry is a function.
#DifEqs calls each function in self.parameterAdjustFuncs if it is not empty, and then each
#of the functions for each differential equation
#This setup allows for an instance of the class to have differential equations swapped out,
#or to add time or other dependencies on the parameters. 
#Check TcellDepletion.py for an example of replaceing differential equations or
#functions for parameters

class ColDetModel:
    #initConditionDict is used in importing paramters to ensure that each initial
    #condition is provided, and in the correct order
    initConditionDict = {'M_R' : 0,
                         'M_Ip' : 1,
                         'M_Id' : 2,
                         'M_A' : 3,
                         'T_0' : 4,
                         'T_1' : 5,
                         'T_2' : 6,
                         'T_80' : 7,
                         'T_8' : 8,
                         'T_c' : 9,
                         'F_a' : 10,
                         'I_y' : 11,
                         'I_12' : 12,
                         'I_10' : 13,
                         'I_4' : 14,
                         'B_I' : 15,
                         'B_E' : 16,
                         'C' : 17,
                         'P_1' : 18,
                         'B_L' : 19}
                         
                         
    #timeToMaturity = 5 #days
    
    def __init__(self, useColEqs=1):
        #self.delM_Iarray = []
        #self.M_Iparray = []
        #self.M_Idarray = []
        #self.timeArray = []
        
        #The two variables below were used when the MMP dynamics were not tracked
        #at the beginning of the model run. With self.colDetModel set to 1, the
        #MMP dynamics are on from the beginning.
        self.colDetFlag = 1
        self.tSwitch = -1
        
        #The below variable determines if the MMP dynamics are included or not
        self.useColEqs = useColEqs
        
        #Parameter from old Kirschner paper. Mostly unused at this point
        self.epsilon = 0.0
        self.parameterAdjustFuncs = {}
        
    
    #Below are the default differential equations used in the model, but they can be replaced with new functions in scripts
    
    
    def dM_Rdt(self,y,t,p):
        #Unpacks the variables from array y
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        #Calculate M_I from the produced and dead M_I
        M_I = M_Ip - M_Id
        
        #Calculate B_T from internal and external bacteria
        B_T = B_I + B_E    
    
        ddt = (p['sr_m']
            +p['a_4a'] * (M_A + p['w_2'] * M_I) 
            +p['sr_4b'] * (F_a / (F_a + p['f_8']*I_10 + p['s_4b'])) 
            -p['k_2']*M_R * (B_E / (B_E + p['c_9']))
            -p['k_3']*M_R * (I_y / (I_y + p['f_1']*I_4 + p['s_1']))  
            *(B_T + p['B']*F_a) / (B_T + p['B']*F_a + p['c_8']) 
            #+p['k_4']*M_A * (I_10 / (I_10 + p['s_8'])) 
            -p['u_MR']*M_R)
            
        return ddt
        
    def dM_Ipdt(self,y,t,p):
        
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        #M_I = M_Ip - M_Id
        
        #B_T = B_I + B_E    
    
        ddt = p['k_2']*M_R * (B_E / (B_E + p['c_9']))
            
        return ddt
        
    def dM_Iddt(self,y,t,p):
        
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        M_I = M_Ip - M_Id
        
        #B_T = B_I + B_E  
        
        if M_I != 0 and B_I != 0:
            
            ddt = (p['k_17']*M_I * (B_I**2 / (B_I**2 + (p['N']*M_I)**2+self.epsilon))  
                    +p['k_14a']*M_I * ((T_c + p['w_3']*T_1)/M_I) / ((T_c + p['w_3']*T_1)/M_I + p['c_4']) 
                    +p['k_14b']*M_I * F_a/(F_a + p['f_9']*I_10 + p['s_4b']) 
                    +p['k_52']*M_I * (( (T_c*( T_1/(T_1+p['c_T1']) ) + p['w_1']*T_1 ) /M_I) 
                    /(( T_c * (T_1/(T_1+p['c_T1']))+p['w_1']*T_1)/M_I + p['c_52'])) 
                    +p['u_MI']*M_I)
            
        else:
            
            ddt = (p['k_17']*M_I * (B_I**2 / (B_I**2 + (p['N']*M_I)**2+self.epsilon + sys.float_info.epsilon))  
                    +p['k_14a']*M_I * ((T_c + p['w_3']*T_1)/(M_I+sys.float_info.epsilon)) / ((T_c + p['w_3']*T_1)/(M_I+sys.float_info.epsilon) + p['c_4']) 
                    +p['k_14b']*M_I * F_a/(F_a + p['f_9']*I_10 + p['s_4b']) 
                    +p['k_52']*M_I * (( (T_c*( T_1/(T_1+p['c_T1']) ) + p['w_1']*T_1 ) /(M_I+sys.float_info.epsilon)) 
                    /(( T_c * (T_1/(T_1+p['c_T1']))+p['w_1']*T_1)/(M_I+sys.float_info.epsilon) + p['c_52'])) 
                    +p['u_MI']*M_I)
            
        return ddt
        
    def dM_Adt(self,y,t,p):
        
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        #M_I = M_Ip - M_Id
        
        B_T = B_I + B_E    
    
        ddt = (p['k_3']*M_R * I_y/(I_y + p['f_1']*I_4 + p['s_1'])
                *(B_T + p['B']*F_a)/(B_T + p['B']*F_a + p['c_8']) 
                -p['k_4']*M_A * (I_10 / (I_10 + p['s_8'])) 
                -p['u_MA']*M_A)
            
        return ddt
        
    def dT_0dt(self,y,t,p):
        
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        M_I = M_Ip - M_Id
        
        #B_T = B_I + B_E    
    
        ddt =  (p['a_1a']*(M_A + p['w_2']*M_I) 
                +p['sr_1b'] * F_a/(F_a + p['f_8'] * I_10 + p['s_4b2'])
                +p['a_2']*T_0 * M_A/(M_A + p['c_15'])
                #In the paper, term 4 is inconsistant with it's corresponding term (term 3 of dT_1/dt) in the paper.
                -p['k_6']*I_12*T_0 * I_y/(I_y + (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
                #Later equations in the paper suggest the form of the term here is correct.
                #-p['k_6']*I_12*T_0 * I_y/(I_y * (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
                -p['k_7']*T_0 * I_4 / (I_4 + p['f_2']*I_y + p['s_2'])
                -p['u_T0']*T_0)
            
        return ddt
        
    def dT_1dt(self,y,t,p):
        
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        M_I = M_Ip - M_Id
        
        #B_T = B_I + B_E    
    
        ddt = (p['a_3a'] * (M_A + p['w_2']*M_I)
                +p['sr_3b'] * F_a/(F_a + p['f_8'] * I_10 + p['s_4b1'])
                #Below is as written in the paper (potential typo)
                +p['k_6']*I_12*T_0 * I_y/(I_y + (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
                #Below agrees with the corresponding term in other equations
                #+p['k_6']*I_12*T_0 * I_y/(I_y*(p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
                -p['u_Ty'] * I_y/(I_y + p['c']) * T_1*M_A 
                -p['u_T1']*T_1)
            
        return ddt
        
    def dT_2dt(self,y,t,p):
        
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        M_I = M_Ip - M_Id
        
        #B_T = B_I + B_E    
    
        ddt = (p['a_3a2'] * (M_A + p['w_2']*M_I)
                +p['sr_3b2'] * F_a/(F_a + p['f_8'] * I_10 + p['s_4b1'])
                +p['k_7']*T_0 * I_4 / (I_4 + p['f_2']*I_y + p['s_2'])
                -p['u_T2']*T_2)
            
        return ddt
        
    def dT_80dt(self,y,t,p):
        
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        M_I = M_Ip - M_Id
        
        #B_T = B_I + B_E    
    
        ddt = (p['a_1a']*(M_A + p['w_2']*M_I)#eqns say wM_I but no w exists, just w1, w2, w3
                +p['sr_1b'] * F_a/(F_a + p['f_8'] * I_10 + p['s_4b2'])
                +p['a_2']*T_80 * M_A/(M_A + p['c_15'])
                #-p['k_6']*I_12*T_80 * I_y/(I_y * (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1']) +
                -p['k_6']*I_12*T_80 * I_y/(I_y + (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
                -p['u_T80']*T_80)
            
        return ddt
        
    def dT_8dt(self,y,t,p):
        
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        M_I = M_Ip - M_Id
        
        #B_T = B_I + B_E    
    
        ddt = (p['m'] * p['a_3ac']*(M_A + p['w_2']*M_I)
                +p['m'] * p['sr_3bc'] * F_a/(F_a + p['f_8'] * I_10 + p['s_4b1'])
                #+p['m'] *p['k_6']*I_12*T_80 * I_y/(I_y * (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
                +p['m'] *p['k_6']*I_12*T_80 * I_y/(I_y + (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
                -p['u_Tcy'] * I_y/(I_y + p['c_c']) * T_8 * M_A
                -p['u_T8']*T_8)
            
        return ddt
        
    def dT_cdt(self,y,t,p):
        
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        M_I = M_Ip - M_Id
        
        #B_T = B_I + B_E    
    
        ddt = (p['m'] * p['a_3ac']*(M_A + p['w_2']*M_I)
                +p['m'] * p['sr_3bc'] * F_a/(F_a + p['f_8'] * I_10 + p['s_4b1'])
                #+p['m'] *p['k_6']*I_12*T_80 * I_y/(I_y * (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
                +p['m'] *p['k_6']*I_12*T_80 * I_y/(I_y + (p['f_1']*I_4 + p['f_7']*I_10) + p['s_1'])
                -p['u_Tcy'] * I_y/(I_y + p['c_c']) * T_c * M_A +
                -p['u_Tc']*T_c)#Steve had u_T8 instead of u_Tc
            
        return ddt
    
    def dF_adt(self,y,t,p):
        
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        M_I = M_Ip - M_Id
        
        B_T = B_I + B_E    
    
        ddt = (p['a_30']*M_I
                +p['a_31']*M_A * (I_y + p['B_2']*B_T)/(I_y + p['B_2']*B_T + (p['f_1']*I_4 + p['f_7']*I_10) + p['s_10'])
                +p['a_32']*T_1
                +p['a_33'] * (T_c + T_8)
                #In the paper, the parameter is written as u_Fa, which isn't included in table III    
                #term5 = -p['u_Fa']*F_a
                #I suspect the subscript is supposed to be u_TNF, as is written below
                - p['u_TNF']*F_a)
            
        return ddt
    
    def dI_ydt(self,y,t,p):
        
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        M_I = M_Ip - M_Id
        
        B_T = B_I + B_E    
    
        ddt = (p['s_g'] * B_T/(B_T + p['c_10']) * I_12/(I_12 + p['s_7']) 
                +p['a_5a']*T_1 * M_A/(M_A + p['c_5a']) 
                +p['a_5b']*T_8 * M_A/(M_A + p['c_5b']) 
                +p['a_5c']*M_I 
                +p['a_7']*T_0 * I_12/(I_12 + p['f_4']*I_10 + p['s_4'])#* M_A/(M_A + 7000.0)
                +p['a_7']*T_80 * I_12/(I_12 + p['f_4']*I_10 + p['s_4'])#* M_A/(M_A + 7000.0) 
                -p['u_iy']*I_y)
            
        return ddt    
    
    def dI_12dt(self,y,t,p):
        
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        #M_I = M_Ip - M_Id
        
        B_T = B_I + B_E    
    
        ddt = (p['s_12'] * B_T/(B_T + p['c_230']) 
                +p['a_23']*M_R * B_T/(B_T + p['c_23'])
                +p['a_8']*M_A * p['s']/(p['s'] + I_10)
                -p['u_i12']*I_12)
            
        return ddt
        
    def dI_10dt(self,y,t,p):
        
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        #M_I = M_Ip - M_Id
        
        #B_T = B_I + B_E    
    
        ddt = (p['d_7']*M_A * p['s_6']/(I_10 + p['f_6']*I_y + p['s_6']) 
                +p['a_16']*T_1
                +p['a_17']*T_2
                +p['a_18']*(T_8 + T_c)
                -p['u_i10']*I_10)
            
        return ddt
        
    def dI_4dt(self,y,t,p):
        
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        #M_I = M_Ip - M_Id
        
        #B_T = B_I + B_E    
    
        ddt = (p['a_11']*T_0
                +p['a_12']*T_2 
                -p['u_i4']*I_4)
            
        return ddt
        
    def dB_Idt(self,y,t,p):
        
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        M_I = M_Ip - M_Id + 1
        
        #B_T = B_I + B_E    

        if M_I != 0 and B_I != 0:
        
            ddt = (p['a_19']*B_I * (1 - B_I**2/(B_I**2 + (p['N']*M_I)**2+self.epsilon))
                +p['k_2'] * p['N']/2 * M_R * B_E/(B_E + p['c_9'])
                -p['k_17']*p['N']*M_I * B_I**2/(B_I**2 + (p['N']*M_I)**2+self.epsilon)
                -p['k_14a']*p['N'] * M_I * ( (T_c + p['w_3']*T_1) / M_I) / 
                ( (T_c + p['w_3']*T_1) / M_I + p['c_4']) 
                -p['k_14b']*p['N']*M_I * F_a / (F_a + p['f_8']*I_10 + p['s_4b'])
                -p['k_52']*p['N']*M_I * (T_c * T_1/(T_1 + p['c_T1']) + p['w_1']*T_1)/M_I /
                ( (T_c * T_1/(T_1 + p['c_T1']) + p['w_1']*T_1 ) /M_I + p['c_52'])
                -p['u_i']*B_I)
            
        else:
            ddt = (p['a_19']*B_I * (1 - B_I**2/(B_I**2 + (p['N']*M_I)**2+self.epsilon + sys.float_info.epsilon))
                +p['k_2'] * p['N']/2 * M_R * B_E/(B_E + p['c_9'])
                -p['k_17']*p['N']*M_I * B_I**2/(B_I**2 + (p['N']*M_I)**2+self.epsilon + sys.float_info.epsilon)
                -p['k_14a']*p['N'] * M_I * ( (T_c + p['w_3']*T_1) / (M_I+sys.float_info.epsilon)) / 
                ( (T_c + p['w_3']*T_1) / (M_I+sys.float_info.epsilon) + p['c_4']) 
                -p['k_14b']*p['N']*M_I * F_a / (F_a + p['f_8']*I_10 + p['s_4b'])
                -p['k_52']*p['N']*M_I * (T_c * T_1/(T_1 + p['c_T1']) + p['w_1']*T_1)/M_I /
                ( (T_c * T_1/(T_1 + p['c_T1']) + p['w_1']*T_1 ) /M_I + p['c_52'])
                -p['u_i']*B_I)
            
        return ddt
        
    def dB_Edt(self,y,t,p):
        
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        M_I = M_Ip - M_Id
        
        #B_T = B_I + B_E    
        
        if C>p['C_L']:
            C = p['C_L']

        if M_I != 0 and B_I != 0:
            
            ddt = (p['a_20']*B_E +
                p['u_i']*B_I +
                -p['k_15']*M_A*B_E +
                -p['k_18']*M_R*B_E +
                p['k_17']*p['N']*M_I * B_I**2 / (B_I**2 + (p['N']*M_I)**2+self.epsilon) +
                -p['k_2'] * p['N']/2 * M_R * B_E/(B_E + p['c_9']) +
                (p['k_14a']*p['N'] * p['N_fracc']*M_I * ((T_c + p['w_3']*T_1) / M_I) / 
                ((T_c + p['w_3']*T_1) / M_I + p['c_4'])) +
                p['k_14b']*p['N'] * p['N_fraca']*M_I * F_a / (F_a + p['f_9']*I_10 + p['s_4b']) +
                -p['k_CL']*p['s_L']*B_E * (1 - C/p['C_L'])/(C + p['s_L']))
                #-p['k_CL']*p['s_L1']*B_E * (p['C_L'] - C)/((p['C_L'] + p['s_L2'])*(C + p['s_L2']))*self.colDetFlag)

            
            
        else:
            ddt = (p['a_20']*B_E +
                p['u_i']*B_I +
                -p['k_15']*M_A*B_E +
                -p['k_18']*M_R*B_E +
                p['k_17']*p['N']*M_I * B_I**2 / (B_I**2 + (p['N']*M_I)**2+self.epsilon+ sys.float_info.epsilon) +
                -p['k_2'] * p['N']/2 * M_R * B_E/(B_E + p['c_9']) +
                (p['k_14a']*p['N'] * p['N_fracc']*M_I * ((T_c + p['w_3']*T_1) / (M_I+sys.float_info.epsilon)) / 
                ((T_c + p['w_3']*T_1) / (M_I+sys.float_info.epsilon) + p['c_4'])) +
                p['k_14b']*p['N'] * p['N_fraca']*M_I * F_a / (F_a + p['f_9']*I_10 + p['s_4b']) +
                -p['k_CL']*p['s_L']*B_E * (1 - C/p['C_L'])/(C + p['s_L']))
                #-p['k_CL']*p['s_L1']*B_E * (p['C_L'] - C)/((p['C_L'] + p['s_L2'])*(C + p['s_L2']))*self.colDetFlag)
            
        return ddt
        
    def dCdt(self,y,t,p):
        
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        #M_I = M_Ip - M_Id
        
        #B_T = B_I + B_E    
    
        ddt = (-p['k_C'] * C/(C + p['k_M']) * P_1 +
                p['sr_C'])
                
        if (C >= p['C_L'] and ddt > 0):
                ddt = 0
            
        return ddt
        
    def dP_1dt(self,y,t,p):
        
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        M_I = M_Ip - M_Id
        
        #B_T = B_I + B_E    
    
        if (self.colDetFlag == 0 or (not self.useColEqs)):
            
            ddt = (-p['u_P1']*P_1 +
                    p['sr_P1'])
 
        else:

            ddt = (-p['u_P1']*P_1 +
                    p['sr_P1'] + 
                    p['B_P1']*M_I+
                    p['a_P1']*M_I * F_a/(F_a + p['s_P1']))
                    #Change M_I to M_Imature
                    #M_I maturity in 5 days?

        return ddt
        
    def dB_Ldt(self,y,t,p):
        
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        M_I = M_Ip - M_Id
        
        #B_T = B_I + B_E    
    
        if (self.colDetFlag == 0 or (not self.useColEqs)):
            
            ddt = 0
            #print(self.colDetFlag)
            
        else:
            
            if C>p['C_L']:
                C = p['C_L']
    
            ddt = p['k_CL']*p['s_L']*B_E * (1 - C/p['C_L'])/(C + p['s_L'])

        return ddt
    
    #DifEqs is the method passed to odeint
    def DifEqs(self,y,t,p):
    
        #Unpacks the variables from array y
        M_R, M_Ip, M_Id, M_A, T_0, T_1, T_2, T_80, T_8, T_c, F_a, I_y, I_12, I_10, I_4, B_I, B_E, C, P_1, B_L = y
        
        #Calculate M_I from the produced and dead M_I
        #M_I = M_Ip - M_Id
        
        #Calculate B_T from internal and external bacteria
        B_T = B_I + B_E
        
        
        
        #if (self.delM_Iarray == []):
        #    self.delM_Iarray.append(M_Ip)
        #else:
        #    self.delM_Iarray.append(M_Ip - self.delM_Iarray[-1])
        #    pass 
        #self.M_Iparray.append(M_Ip)
        #self.M_Idarray.append(M_Id)
        #self.timeArray.append(t)
        
        #5 days to maturity
        try:
            if self.parameterAdjustFuncs != {} and self.baseParameters:
                for key, value in self.parameterAdjustFuncs.items():
                    p[key] = value(y,t,p, self.baseParameters)
        except:
            pass
        
        dM_Rdt = self.dM_Rdt(y,t,p)
                
        dM_Ipdt = self.dM_Ipdt(y,t,p)
                
        dM_Iddt = self.dM_Iddt(y,t,p)
                
        dM_Adt = self.dM_Adt(y,t,p)
        
        dT_0dt = self.dT_0dt(y,t,p)
                
        dT_1dt = self.dT_1dt(y,t,p)
                
        dT_2dt = self.dT_2dt(y,t,p)
    
        dT_80dt = self.dT_80dt(y,t,p)
    
        dT_8dt = self.dT_8dt(y,t,p)
                
        dT_cdt = self.dT_cdt(y,t,p)
                
        dF_adt = self.dF_adt(y,t,p)
                
        dI_ydt = self.dI_ydt(y,t,p)
                   
        dI_12dt = self.dI_12dt(y,t,p)
                
        dI_10dt = self.dI_10dt(y,t,p)
                
        dI_4dt = self.dI_4dt(y,t,p)
                
        dB_Idt = self.dB_Idt(y,t,p)
                
        dB_Edt = self.dB_Edt(y,t,p)
                
        dCdt = self.dCdt(y,t,p)
        
        dP_1dt = self.dP_1dt(y,t,p)
        
        dB_Ldt = self.dB_Ldt(y,t,p)
        
        #The below was used to switch on the MMP dynamic when a steady state was reached            
        try:
            dydtArray = np.array([dM_Rdt, dM_Ipdt, dM_Iddt, dM_Adt, dT_0dt, dT_1dt, dT_2dt, dT_80dt, dT_8dt, dT_cdt, dF_adt, dI_ydt, dI_12dt, dI_10dt, dI_4dt, dB_Idt, dB_Edt, dCdt, dP_1dt])
            dydt_y = np.append([dydtArray[0]/y[0], (dydtArray[1]-dydtArray[2])/(y[1]-y[2])], dydtArray[3:19]/y[3:19])
            
            #Switches on when every derivative is %0.002 of the respective variable
            if (np.all(np.less_equal(dydt_y, 0.00002)) and (not self.colDetFlag) and self.useColEqs):
                self.colDetFlag = 1
                self.tSwitch = t
        except ZeroDivisionError:
            pass
        
        return [dM_Rdt, dM_Ipdt, dM_Iddt, dM_Adt, dT_0dt, dT_1dt, dT_2dt, dT_80dt, dT_8dt, dT_cdt, dF_adt, dI_ydt, dI_12dt, dI_10dt, dI_4dt, dB_Idt, dB_Edt, dCdt, dP_1dt,dB_Ldt]

           