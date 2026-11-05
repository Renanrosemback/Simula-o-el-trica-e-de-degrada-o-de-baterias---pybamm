import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pybamm


from battery_circuit_POO import PerfilPotencia

import config
import math
import time
import winsound




#========================================================================================
#POTÊNCIA
#========================================================================================
segundos_simulados = 2*1590 #5horas 

pot = PerfilPotencia()
#P_unico = pot.leitura_csv() #entra no STROE, porque é modulado lá dentro
P_unico = pot.carregar_perfil_potencia(segundos_totais=segundos_simulados)
P_uni = -P_unico/(config.Nsb*1000) #entra no pybamm porque deve ser modulado antes


#pybamm:
#positivo - out
#negativo - in

#stroe:
#positivo - in
#negativo - out


class Bateria:
    DT = config.DT  # s - passo fixo, igual ao FixedStep=1 do .plecs original
    OCV_TABLE = config.OCV_TABLE #valores LUT no PLECS
    SOC_STEP = config.SOC_STEP #medidos a partir de 5% de SOC = 0.05

    #alterado:
    SOCmax=80
    SOCmin=20
    SOCi0 = config.SOCi0
    CAh0=config.CAh0
    Cr = config.Cr
    Nsb=config.Nsb

    def __init__(self, Rs=config.Rs0, R1=config.R10,
               R2=config.R20, Q1=config.Q10, Q2=config.Q20, 
               Npb=config.Npb_padrao):
        self.Rs = Rs 
        self.R1 = R1 
        self.R2 = R2 
        self.Q1 = Q1
        self.Q2 = Q2 
        self.CAh = Bateria.CAh0
        self.soc = Bateria.SOCi0 / 100.0 #valor do soc de 0 até 1
        
        self.soc_max_frac = Bateria.SOCmax / 100.0 #0.8
        self.soc_min_frac = Bateria.SOCmin / 100.0 #0.2
        self.Ilim = Bateria.CAh0 * Bateria.Cr 
        
        self.Npb = Npb
        self.Nsb = Bateria.Nsb #não pode ser alterada
    
    def _ocv(self): 
        """
        A partir da não linearidade da bateria faz uma interpolação para 
        obter o valor que será utilizado como tensão de circuito aberto. """
        SOC_STEPS = len(Bateria.OCV_TABLE) - 1  #20 
        idx_f = self.soc / Bateria.SOC_STEP
        if idx_f <= 0.0:
            return Bateria.OCV_TABLE[0]
        elif idx_f >= SOC_STEPS:
            return Bateria.OCV_TABLE[SOC_STEPS]
        idx0 = int(idx_f)
        frac = idx_f - idx0
        return Bateria.OCV_TABLE[idx0] + frac * (Bateria.OCV_TABLE[idx0 + 1]
                                                 - Bateria.OCV_TABLE[idx0])
    
    def calc_zarc_params(self):
        """
        Recalcula os parametros dos dois ramos ZARC (Ca/Ra/Cb/Rb/Cc/Rc/Cd/Rd/Ce/Re,
        para os ramos 1 e 2) a partir dos valores correntes de R1, Q1, R2, Q2.
        Espelha exatamente as contas que estao tanto no SOCestimation.m (dentro
        do loop, apos atualizar R1/R2) quanto no calculo inicial do SOCestimation.py.
        """
        w_zarc1 = (1 / (self.Q1 * self.R1)) ** (1 / 0.5)
        Rc1 = 0.906 *self.R1 * (np.sin(0.5 * np.pi * 0.5 / 2)) / (1 + np.cos(0.5 * np.pi * 0.5 / 2))
        Cc1 = 1 / (w_zarc1 * Rc1)
        Rb1, Rd1 = 0.670 * (self.R1 - Rc1) / 2, 0.670 * (self.R1 - Rc1) / 2
        Cb1 = 1 / (w_zarc1 * 10.886 * Rb1)
        Cd1 = 1 / (w_zarc1 * (1 / 10.886) * Rd1)
        Ra1, Re1 = (self.R1 - Rc1 - Rb1 - Rd1) / 2, (self.R1 - Rc1 - Rb1 - Rd1) / 2
        Ca1 = 1 / (w_zarc1 * (0.67 ** 2) * Ra1)
        Ce1 = 1 / (w_zarc1 * (1 / 0.67 ** 2) * Re1)

        w_zarc2 = (1 / (self.Q2 * self.R2)) ** (1 / 0.5)
        Rc2 = 0.833 * self.R2 * (np.sin(0.8 * np.pi * 0.5 / 2)) / (1 + np.cos(0.8 * np.pi * 0.5 / 2))
        Cc2 = 1 / (w_zarc2 * Rc2)
        Rb2, Rd2 = 0.804 * (self.R2 - Rc2) / 2, 0.804 * (self.R2 - Rc2) / 2
        Cb2 = 1 / (w_zarc2 * 4.675 * Rb2)
        Cd2 = 1 / (w_zarc2 * (1 / 4.675) * Rd2)
        Ra2, Re2 = (self.R2 - Rc2 - Rb2 - Rd2) / 2, (self.R2 - Rc2 - Rb2 - Rd2) / 2
        Ca2 = 1 / (w_zarc2 * (0.804 ** 2) * Ra2)
        Ce2 = 1 / (w_zarc2 * (1 / 0.804 ** 2) * Re2)
        zarc = dict(Ca1=Ca1, Ra1=Ra1, Cb1=Cb1, Rb1=Rb1, Cc1=Cc1, Rc1=Rc1, Cd1=Cd1, Rd1=Rd1, Ce1=Ce1, Re1=Re1,
                    Ca2=Ca2, Ra2=Ra2, Cb2=Cb2, Rb2=Rb2, Cc2=Cc2, Rc2=Rc2, Cd2=Cd2, Rd2=Rd2, Ce2=Ce2, Re2=Re2)
        return zarc
    
    def _zoh_coeffs(self, R, C, dt=config.DT):
        """
        Coeficientes de discretizacao exata (zero-order hold, corrente
        constante dentro do passo) de um ramo RC paralelo: v[n+1] = a*v[n] + b*I[n]."""
        tau = R * C
        a = math.exp(-dt / tau)
        b = R * (1.0 - a)
        return a, b

    def _branch_coeffs(self): 
        """
        Monta as listas de coeficientes (a, b) para os 5 sub-ramos de ZARC1 e
        ZARC2, na ordem Ra/Ca, Rb/Cb, Rc/Cc, Rd/Cd, Re/Ce (mesma ordem em serie
        do schematic ZARC1/ZARC2 do .plecs)."""
        zarc = self.calc_zarc_params()
        ramos1 = [(zarc['Ra1'], zarc['Ca1']), (zarc['Rb1'], zarc['Cb1']),
                  (zarc['Rc1'], zarc['Cc1']), (zarc['Rd1'], zarc['Cd1']),
                  (zarc['Re1'], zarc['Ce1'])]
        ramos2 = [(zarc['Ra2'], zarc['Ca2']), (zarc['Rb2'], zarc['Cb2']),
                  (zarc['Rc2'], zarc['Cc2']), (zarc['Rd2'], zarc['Cd2']),
                  (zarc['Re2'], zarc['Ce2'])]
        coefs1 = [self._zoh_coeffs(R, C) for R, C in ramos1] #joga cada parte da tupla 
        coefs2 = [self._zoh_coeffs(R, C) for R, C in ramos2]
        
        return coefs1, coefs2

    def _calcula_corrente(self, pot_t, V_bess_prev):
        """
        Define a corrente a partir da tensao de acordo com os limites
        de corrente nominais de saturacao. Além de lmitar o máximo de corrente que 
        pode ser retirada devido a saturação pelo soc"""
        #alteracao: retirei o Npb
        I = pot_t / (V_bess_prev) 
        
        if I > self.Ilim:
            I = self.Ilim
        elif I < -self.Ilim:
                I = -self.Ilim

        I_min_soc = ((self.soc_min_frac - self.soc)*self.CAh*3600/self.DT)

        I_max_soc = ((self.soc_max_frac - self.soc)*self.CAh * 3600 / self.DT)

        # Corrente efetivamente aplicada
        I = np.clip(I, I_min_soc, I_max_soc)

        return I

    def _atualiza_tensao_ramos_zarc(self, v1, v2, coefs1, coefs2, I):
        for k in range(5):
            a, b = coefs1[k]
            v1[k] = a * v1[k] + b * I
    
            a, b = coefs2[k]
            v2[k] = a * v2[k] + b * I
            
        
        return v1, v2

    def _calcula_tensao(self, ocv, I, sum_v):
            """
            Serve para obter a tensão Vbees a partir do zarc, fonte não ideal ocv e 
            resistencia Rs
            """
            v_bess = ocv + I * self.Rs + sum_v
            return v_bess


    def simular_celula(self, P_pack):
        """
        Simula a dinamica eletrica (Rs + ZARC1 + ZARC2) e a integracao de SOC de
        uma string de celulas ao longo de `len(P_pack)` segundos, com os ramos
        RC comecando relaxados (v=0) a cada chamada - mesma condicao inicial que
        o .plecs original usava a cada `plecs.simulate()` por mes (v_init=0 nos
        capacitores; nenhum estado de relaxacao eletrica era carregado de um mes
        para o outro, so R/C/CAh/SOCi, via os parametros do modelo).

        P_pack: sequencia de potencia total do pack (W), convencao positivo =
            carga/regenerativa, negativo = descarga/tracao (mesmo sinal de
            mission_day_pot em carregar_perfil_potencia).
        zarc: dict com Ca1..Re1, Ca2..Re2 (saida de calc_zarc_params).
        CAh0: capacidade nominal de fabrica (Ah) - define o limite de corrente
            (Cr-rate sobre a capacidade nominal, nao a degradada).
        CAh: capacidade atual, ja degradada (Ah) - usada na integracao do SOC.
        SOCmax, SOCmin: limites de SOC em % (0-100).
        SOCi0: SOC inicial em % (0-100).

        Retorna (SOC_pct, P_bess_filtrada): arrays (listas) de tamanho
        len(P_pack). SOC_pct em % (0-100), amostrado no INICIO de cada passo
        (SOC_pct[0] == SOCi0). P_bess_filtrada em W, potencia total do pack
        (ja multiplicada por Npb), zerada quando o SOC de inicio de passo esta
        fora do intervalo aberto (SOCmin, SOCmax).
        """
        #alteração: P_unico fica mais simples, já que imaginamos
        #que a potência dividida para cada célula é a mesma, fica mais 
        #fácil já deixar isso explicito.
        P_unico = P_pack/(self.Npb*self.Nsb)
        n = len(P_unico)
        coefs1, coefs2 = self._branch_coeffs()
        
        v1 = [0.0, 0.0, 0.0, 0.0, 0.0]
        v2 = [0.0, 0.0, 0.0, 0.0, 0.0]

        #Definicao variaveis de saida:
        SOC_pct = [0.0] * n
        P_bess_filtrada = [0.0] * n
        I_hist = [0.0] * n 
        V_hist = [0.0] * n 
        
        ocv = self._ocv()
        v_bess = ocv 
        sum_v = sum(v1) + sum(v2)

        #condicoes iniciais: #supondo antes da ativacao da bateria
        SOC_pct[0] = self.soc * 100
        P_bess_filtrada[0] = 0
        I_hist[0] = 0
        V_hist[0] = v_bess #v_bess[i-1] 
        
        for i in range(1, n):
            pot_t = P_unico[i] #variavel de entrada

            I = self._calcula_corrente(pot_t, v_bess) #I[i] - a partir de v_bess[i-1]
            
            self.soc += I / (self.CAh * 3600.0) #soc[i]
            
            v1, v2 = self._atualiza_tensao_ramos_zarc(v1, v2, coefs1, coefs2, I) #zarc[i]

            ocv = self._ocv() 
            sum_v = sum(v1) + sum(v2) #ovc[i]

            v_bess = self._calcula_tensao(ocv, I, sum_v) #v_bess[i]

            P_bess = I * self.Npb * self.Nsb * v_bess #p_bess[i]

            #salvando resultados
            I_hist[i] = I
            V_hist[i] = v_bess
            SOC_pct[i] = self.soc * 100.0
            P_bess_filtrada[i] = P_bess
                
        return SOC_pct, P_bess_filtrada, I_hist, V_hist
    
    def aplicar_degradacao(
        self,
        Rs_i,
        R1_i,
        R2_i,
        Q1_i,
        Q2_i,
        C_fade_cyc,
        C_fade_cal, 
        **kwargs
    ):
        self.Rs += Rs_i
        self.R1 += R1_i
        self.R2 += R2_i
        self.Q1 -= Q1_i
        self.Q2 -= Q2_i
        self.CAh = Bateria.CAh0 * (1 - (C_fade_cyc + C_fade_cal)/100)
        
    def resetar(self):
        self.Rs = config.Rs0
        self.R1 = config.R10
        self.R2 = config.R20
        self.Q1 = config.Q10
        self.Q2 = config.Q20
    
        self.CAh = Bateria.CAh0
        self.soc = Bateria.SOCi0 / 100.0
    

#========================================================================================
#CRIAÇÃO DE UMA CLASSE 
#========================================================================================

class BateriaPyBaMM():
    Nsb = config.Nsb #configurado em 390 células - tensão nominal de 1287 do pack
    soc_inicial = config.SOCi0/100
    SOCmin = config.SOCmin
    SOCmax = config.SOCmax

    def __init__(self, celula='A123 26650M1', modelo='SPM', degradacao=False):
        print(pybamm.__version__) #a versão é importante: 26.9.0.0

        self.Parametro_celula_correta(celula, degradacao)

        self.Nsb = BateriaPyBaMM.Nsb
        self.soc_inicial = BateriaPyBaMM.soc_inicial

        self.SOCmin = BateriaPyBaMM.SOCmin/100
        self.SOCmax = BateriaPyBaMM.SOCmax/100

        
                        
        self.CAh0 = float(self.param["Nominal cell capacity [A.h]"])
        self.Ilim = self.CAh0 * config.Cr
        


    def Parametro_celula_correta(self, celula, degradacao):

        self.param = pybamm.ParameterValues("Prada2013")
        if celula == 'A123 26650M1B':
            ocp_original = self.param["Positive electrode OCP [V]"]
            j0_original = self.param["Negative electrode exchange-current density [A.m-2]"]
            kappa_original = self.param["Electrolyte conductivity [S.m-1]"]
            D_original = self.param["Electrolyte diffusivity [m2.s-1]"]
            def ocp_a123(*args):
                return ocp_original(*args) + 0.02515325

            def j0_a123(*args):
                return 3.53854 * j0_original(*args)

            def kappa_a123(*args):
                return 0.20 * kappa_original(*args)

            D_a123 = 1.06384 * D_original

            self.param.update({
                    "Positive electrode OCP [V]": ocp_a123,
                    "Negative electrode exchange-current density [A.m-2]": j0_a123,
                    "Electrolyte conductivity [S.m-1]": kappa_a123,
                    "Electrolyte diffusivity [m2.s-1]": D_a123,})
        if degradacao:
            self.param.update({
                    "Initial SEI thickness [m]": 5e-9,
                    "SEI solvent diffusivity [m2.s-1]": 8.84e-20,
                    "Bulk solvent concentration [mol.m-3]": 4541.0,
                    "SEI partial molar volume [m3.mol-1]": 0.162 / 1690,
                    "SEI growth activation energy [J.mol-1]": 55500.0,
                    "SEI resistivity [Ohm.m]": 1 / (17.5e-5),
                    "Ratio of lithium moles to SEI moles": 2.0,  # Hipótese estequiométrica
                        }, check_already_exists=False)

  
    def resultados_SPM(self):

        corrente = self.solution["Current [A]"].entries
        tensao = self.solution["Voltage [V]"].entries
        tempop = self.solution["Time [s]"].entries
        perda_CAh = self.solution["Loss of capacity to negative SEI [A.h]"].entries

        CAh_atual = self.CAh0 - perda_CAh
        
        return corrente, tensao, tempop, perda_CAh, CAh_atual

    def _experiment_perfil_pot(self, P_uni):
        """O objetivo desse método é simular a célula a partir do perfil de potência e 
        criar uma função limitadora de corrente. """
        
        faixa_soc = 0.002
        self.CAh = float(self.param["Nominal cell capacity [A.h]"])

        t = np.arange(0, len(P_uni))
        P_uni = np.column_stack((t, P_uni))
        
        potencia = pybamm.Interpolant(t, P_uni[:, 1], pybamm.t, interpolator="linear")

        def controle_corrente(variables):
            I = variables["Current [A]"]
            V = variables["Voltage [V]"]

            # Calculando SOC
            Q_desc = variables["Discharge capacity [A.h]"]
            soc = self.soc_inicial - Q_desc / self.CAh

            #Saturando corrente
            I_solicitada = potencia / V
            I_limitada = pybamm.minimum(pybamm.maximum(I_solicitada, -self.Ilim),self.Ilim)

            #Saturando soc em 20%
            fator_descarga = pybamm.minimum(1,pybamm.maximum(0, (soc - self.SOCmin) / faixa_soc))

            #Saturando soc em 80%
            fator_carga = pybamm.minimum(1, pybamm.maximum(0, (self.SOCmax - soc) / faixa_soc))

            I_descarga = pybamm.maximum(I_limitada, 0)
            I_carga = pybamm.minimum(I_limitada, 0)

            # NOVO: Corrente efetivamente aplicada
            I_aplicada = (I_descarga * fator_descarga + I_carga * fator_carga)

            return I - I_aplicada

        experiment = pybamm.Experiment([
                pybamm.step.CustomStepImplicit(controle_corrente, duration=float(t[-1]), period=1)])
        return experiment
        

    def simulacao_SPM(self, P_uni):
        inicio = time.perf_counter()
        

        options = {
            "SEI": "solvent-diffusion limited",
            "SEI film resistance": "average",
        }

        experiment = self._experiment_perfil_pot(P_uni)
        
        model = pybamm.lithium_ion.SPM(options=options)#options=options

        self.sim = pybamm.Simulation(model=model, experiment=experiment, parameter_values=self.param) #, parameter_values=self.param
        
        self.solution = self.sim.solve(initial_soc=self.soc_inicial)
        
        tempo_execucao = time.perf_counter() - inicio
        print(f'O tempo de execução para uma célula SPM foi: {tempo_execucao}')
        

#========================================================================================
#SIMULAÇÃO COM EEC DO FABRICIO
#========================================================================================

inicio = time.perf_counter() 
pack1 = Bateria(Npb=1000)  #essa classe de bateria está alterada para obter somente o SOC
SOC_pct, P_bess_filtrada, I_hist, V_hist = pack1.simular_celula(P_pack= P_unico)
tempo_execucao = time.perf_counter() - inicio
print(f'O tempo de execução para uma célula do fabrício foi: {tempo_execucao}')

#========================================================================================
#SIMULAÇÃO COM PYBAMM:
#========================================================================================

bateria2 = BateriaPyBaMM(celula='A123 26650M1B', degradacao=True)
bateria2.simulacao_SPM(P_uni)


#========================================================================================
#EXTRAINDO RESULTADOS DAs SIMULAÇÕES:
#========================================================================================

corrente, tensao, tempop, perda_CAh, CAh_atual = bateria2.resultados_SPM()

tempof = np.arange(0, len(P_unico), 1)
V_hist = np.array(V_hist)
I_hist = -np.array(I_hist) #tem que inverter a corrente para ajustar as convenções
P_bess_filtrada = -np.array(P_bess_filtrada) #tem que inverter a corrente para ajustar as convenções

#POTENCIAS CONFERIDAS: 
Pot_stroe = V_hist*I_hist
Pot_spm = tensao*corrente

#ERROS DE CORRENTE:
corrente_ajus = np.interp(tempof, tempop, corrente)
erro_I = 100*(np.abs(I_hist - corrente_ajus)) / I_hist

#ERROS DE TENSÃO:
tensao_ajus = np.interp(tempof, tempop, tensao)
erro_V = 100*(np.abs(V_hist - tensao_ajus)) / V_hist

#EXTRAINDO SOC:
carga_acumulada = np.zeros(len(tempop))

for k in range(1, len(tempop)):
    dt = tempop[k] - tempop[k - 1]
    corrente_media = (corrente[k] + corrente[k - 1]) / 2
    carga_acumulada[k] = (carga_acumulada[k - 1] + corrente_media * dt)

SOC = 100*(bateria2.soc_inicial - carga_acumulada / (3600 * config.CAh0))

winsound.Beep(1000, 500)




#========================================================================================
#PLOTANDO VALORES:
#========================================================================================

#PRINT DA TENSÃO
plt.figure()
plt.plot(tempop, tensao, color='red', label='Tensão SPM')
plt.plot(tempof, V_hist, color='blue',  label='Tensão STROE \nConvenção: BrakingPower negativo')
#plt.plot(tempof, V_hist1, color="#F56E00",  label='Tensão bess ')
#plt.plot(tempof, V_hist2, color='#FFD700',  label='Tensão sum_v')
#plt.plot(tempof, V_hist3, color="#C5A206",  label='Tensão ocv')
plt.xlabel("Tempo [s]")
plt.ylabel("Tensão [V]")
plt.title('Comparação de tensão sem alteração de CH0')
plt.legend(fontsize="small")
plt.grid()
plt.show()

#PRINT DO SOC 
plt.figure()
plt.xlabel("Tempo [s]")
plt.ylabel("SOC")
plt.title("Estado de carga obtido por integração de coulomb")
plt.plot(tempop, SOC, color='red', label='SOC SPM')
plt.plot(tempof, SOC_pct, color='blue', label='SOC STROE')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


#PRINT DA CORRENTE:
plt.figure()
plt.plot(tempop, corrente, color='red', label='Corrente SPM')
plt.plot(tempof, I_hist, color='blue', label='Corrente STROE  \nConvenção: BrakingPower negativo')
plt.xlabel("Tempo [s]")
plt.ylabel("Corrente [A]")
plt.legend(fontsize="small")
plt.title("Comparação de correntes sem alteração de CH0")
plt.grid()
plt.show()


# PLOT DA DEGRADAÇÃO DA CAPACIDADE
plt.figure()
plt.plot(tempop, perda_CAh, color='red', label='Perda de capacidade SPM')
plt.plot(tempop, CAh_atual, color='blue', label='Capacidade atual SPM')
plt.xlabel("Tempo [s]")
plt.ylabel("Capacidade [Ah]")
plt.title('Evolução da capacidade da célula - SPM')
plt.legend(fontsize="small")
plt.grid()
plt.show()

#PLOTANDO EM CONJUNTO VALORES OBTIDOS PELO FABRICIO 


fig, axs = plt.subplots(4, 1, figsize=(12, 10), sharex=True)
# 1. SOC
axs[0].plot(tempof, SOC_pct, color='green')
axs[0].set_ylabel('SOC [%]')
axs[0].set_title('Estado de carga (SOC)')
axs[0].grid()

# 2. Potência
axs[1].plot(tempof, P_bess_filtrada, color='blue')
axs[1].set_ylabel('Potência [W]')
axs[1].set_title('Potência da bateria')
axs[1].grid()

# 3. Corrente
axs[2].plot(tempof, I_hist, color='red')
axs[2].set_ylabel('Corrente [A]')
axs[2].set_title('Corrente da célula')
axs[2].grid()

# 4. Tensão
axs[3].plot(tempof, V_hist, color='orange')
axs[3].set_ylabel('Tensão [V]')
axs[3].set_xlabel('Tempo [s]')
axs[3].set_title('Tensão da célula')
axs[3].grid()

plt.tight_layout()
plt.show()


#OBTENDO ERRO ABSOLUTO DAS CORRENTES:
'''plt.figure()
plt.plot(tempof, erro_I, color='green', label='(I_spm- I_stroe)/I_stroe %')
plt.plot(tempof, I_hist, color='red', label='corrente STROE')
plt.title('Erro absoluto de I dos modelos')
plt.xlabel("Tempo [s]")
plt.ylabel("Corrente [A] e %")
plt.legend(fontsize="small")
plt.grid()
plt.show()'''

#OBTENDO ERRO ABSOLUTO DAS TENSÕES:
'''plt.figure()
plt.plot(tempof, erro_V, color='green', label='(V_spm - V_stroe)/V_stroe %')
plt.plot(tempof, V_hist, color='red', label='tensão STROE')
plt.title('Erro absoluto de V dos modelos')
plt.xlabel("Tempo [s]")
plt.ylabel("Tensão [V] e %")
plt.legend(fontsize="small")
plt.grid()
plt.show()
'''

#grafico erro relativo das correntes
'''erro_I_relativo = np.where(I_stroe != 0, corrente / I_stroe, 0.0)
plt.figure()
plt.plot(tempop, erro_I_relativo, color='yellow', label='(I_SPM - I_STROE)/ I_STROE')
plt.title('Erro relativo em I dos modelos')
plt.xlabel("Tempo [s]")
plt.ylabel("Corrente [A]")
plt.legend(fontsize="small")
plt.grid()
plt.show()'''


#PLOTANDO POTÊNCIA: todos os metodos chegaram na potencia 
'''plt.figure(figsize=(12, 5))
plt.plot(tempof, P_uni, color='red', label='Perfil potência')
plt.plot(tempof, Pot_stroe, color="#F56200", label='I*V STROE')
plt.plot(tempop, Pot_spm, color='blue', label= 'I*V SPM')
plt.xlabel("Tempo [s]")
plt.ylabel("Potência [kW]")
plt.title("Perfil de Potência aplicado ao SPM")
plt.legend(fontsize="small")
plt.grid()
plt.show()

'''

#obtendo o soc a partir da estequiometria:
'''theta_n = solution["X-averaged negative particle stoichiometry"].entries
theta_p = solution["X-averaged positive particle stoichiometry"].entries
theta_n_media = np.mean(theta_n, axis=0) #precisa da media pq ele faz estequiometria em 20 pontos
theta_p_media = np.mean(theta_p, axis=0)
        '''

'''"X-averaged negative particle stoichiometry",
"X-averaged positive particle stoichiometry",'''








