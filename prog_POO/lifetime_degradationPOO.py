import numpy as np
from scipy.signal import find_peaks
import rainflow
import config

from battery_circuit_POO import Bateria
from battery_circuit_POO import PerfilPotencia
import time



class BatteryDegradation:
    Rs0 = config.Rs0
    R10 = config.R10
    R20 = config.R20
    Q10 = config.Q10
    Q20 = config.Q20
    
    def __init__(self, Tbat, C_fade_cyc_in, C_fade_Cal_in):
        
        self.Tbat = Tbat
        self._atualizar_acumulado(C_fade_cyc_in, C_fade_Cal_in)
        
        
    def _atualizar_acumulado(self, C_fade_cyc, C_fade_cal):
        self.C_fade_cyc_in = C_fade_cyc
        self.C_fade_Cal_in = C_fade_cal
        
    def resetar(self):

        self.C_fade_cyc_in = 0.0
        self.C_fade_Cal_in = 0.0
        
    def envelhecimento_calendarico(self, mission_day_SOC):
        diffs = np.diff(mission_day_SOC) #pega a diferença entre vaores. 
        mudou = np.where(diffs != 0)[0] #retorna os pontos que houveram essas variações; 
        inicio_bloco = np.concatenate(([0], mudou + 1))
        fim_bloco = np.concatenate((mudou, [len(mission_day_SOC) - 1]))
        comprimento_bloco = fim_bloco - inicio_bloco + 1
        time_ca1 = (comprimento_bloco[:-1] - 1).astype(float)
        SOC_cal = mission_day_SOC[inicio_bloco[:-1]] #é o vaor de soc que vai para a equação de degradação. 
        tt = np.sum(time_ca1)

        Tbat_mean_wo = self.Tbat
        
        # Cálculo do Capacity Fade de Calendário Mensal
        C_fade_cal_sum = (1.9775e-11) * np.exp(0.07511 * (Tbat_mean_wo + 273.15)) * 1.639 * np.exp(0.007388 * SOC_cal) * ((time_ca1 / (24*30*60*60))**0.8)
                         
        # Acúmulo do Fade de Calendário
        p_cal = 10 / 8
        C_fade_Cal_out = (self.C_fade_Cal_in**p_cal + np.sum(C_fade_cal_sum**p_cal))**(1 / p_cal)

        # Aumentos de Resistência e Alterações no Circuito RC
        #de onde vieram essas equações??
        R_fade_c = BatteryDegradation.Rs0 * (1/100) * 1.6487e-11 * np.exp(0.07559 * (Tbat_mean_wo + 273.15)) * (0.5012 * SOC_cal**(0.4259)) * (time_ca1 / (24*30*60*60))
        Rsf = np.sum(R_fade_c)
        
        R1_fade_c = BatteryDegradation.R10 * (1/100) * 2.1165e-12 * np.exp(0.08192 * (Tbat_mean_wo + 273.15)) * 0.9778 * SOC_cal**(0.4118) * (time_ca1 / (24*30*60*60))
        R1_i = np.sum(R1_fade_c)
        
        Q1_fade_c = BatteryDegradation.Q10 * (1/100) * 7.6043e-11 * np.exp(0.03592 * (Tbat_mean_wo + 273.15)) * 1.673 * SOC_cal**(0.1121) * (time_ca1 / (24*30*60*60))
        Q1_i = np.sum(Q1_fade_c)
        
        return  C_fade_Cal_out, R1_i,  Q1_i, Rsf, tt
    
    def _rainflow(self, mission_day_SOC):
        pos_loc, _ = find_peaks(mission_day_SOC) #pega a posição dos máximos
        neg_loc, _ = find_peaks(-mission_day_SOC) #pega a posição dos mínimos
        
        #pega a vaor dos máximos e mínimos:
        pos_pks = mission_day_SOC[pos_loc]
        neg_pks = mission_day_SOC[neg_loc]
        
        time_locs = np.concatenate((pos_loc, neg_loc))
        peaks_total = np.concatenate((pos_pks, neg_pks))
        
        # Reconstruindo o sinal apenas com picos e vales ordenados no tempo
        sort_idx = np.argsort(time_locs)
        mission_oneyear_SOCC = peaks_total[sort_idx]
        
        # Algoritmo de Rainflow para contar os ciclos
        SOC_cycling = []
        DOD_cycling = []
        
        # rainflow.extract_cycles devolve (rng, mean, count, i_start, i_end)
        for rng, mean, count, i_start, i_end in rainflow.extract_cycles(mission_oneyear_SOCC):
            SOC_cycling.append(mean)
            DOD_cycling.append(rng) 
            
        SOC_cycling = np.array(SOC_cycling)
        DOD_cycling = np.array(DOD_cycling)
        
        return SOC_cycling, DOD_cycling
        
    def envelhecimento_ciclico(self, mission_day_SOC):
        Tb_cycling = self.Tbat
        SOC_cycling, DOD_cycling = self._rainflow(mission_day_SOC)
        
        C_fade_cyc1 = 2.6418 * np.exp(-0.01943 * SOC_cycling) * (0.004 * np.exp(0.01705 * (Tb_cycling + 273.15))) *  (0.0123 * DOD_cycling**(0.7162))
        C_fade_cyc_out = np.sqrt(self.C_fade_cyc_in**2 + np.sum(C_fade_cyc1**2))

        # Aumentos de Resistência e Alterações no Circuito RC por Ciclagem
        Rs_sum_i = BatteryDegradation.Rs0 * (1/100) * 7.799e-14 * np.exp(0.0934 * (Tb_cycling + 273.15)) * 2.356e-5 * DOD_cycling**(0.9395)
        Rsc = np.sum(Rs_sum_i)
        
        
        Rs2_sum_i = BatteryDegradation.R20 * (1/100) * 1.3636e-11 * np.exp(0.07742 * (Tb_cycling + 273.15)) * 5.342e-4 * DOD_cycling
        R2_i = np.sum(Rs2_sum_i)
        
        Qs2_sum_i = BatteryDegradation.Q20 * (1/100) * 9.457e-11 * np.exp(0.0358 * (Tb_cycling + 273.15)) * 3.571e-5 * DOD_cycling
        Q2_i = np.sum(Qs2_sum_i)
        
        return C_fade_cyc_out, Rsc, R2_i, Q2_i

    def lifetime_degradation(self, mission_day_SOC):
        
        mission_day_SOC = np.asarray(mission_day_SOC)
        C_fade_Cal_out, R1_i,  Q1_i, Rsf, tt = self.envelhecimento_calendarico(mission_day_SOC)
        C_fade_cyc_out, Rsc, R2_i, Q2_i = self.envelhecimento_ciclico(mission_day_SOC)
        Rs_i = Rsf + Rsc
        
        return {'Rs_i': Rs_i,'R1_i': R1_i,'R2_i': R2_i, 
                'Q1_i': Q1_i, 'Q2_i': Q2_i,'C_fade_cyc': C_fade_cyc_out,
                'C_fade_cal': C_fade_Cal_out, 'tt': tt}
        



#==========================TESTE DE TEMPO DE DEGRADAÇÃO:=================================

pot = PerfilPotencia()

P_pack = pot.carregar_perfil_potencia()
#Potencia 

inicio = time.perf_counter() 
pack1 = Bateria( Npb=1000) #limita a parada para nao 
SOC_pct, P_bess_filtrada = pack1.simular_celula(P_pack=P_pack)
tempo_execucao = time.perf_counter() - inicio


print(f'O tempo de execução para 1000 células {tempo_execucao}')

inicio = time.perf_counter() 
pack1 = Bateria( Npb=3000) #limita a parada para nao 
SOC_pct, P_bess_filtrada = pack1.simular_celula(P_pack=P_pack)
tempo_execucao = time.perf_counter() - inicio

print(f'O tempo de execução para 3000 células {tempo_execucao}')