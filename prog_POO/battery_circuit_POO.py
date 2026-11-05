import config
import math
import numpy as np
import pandas as pd



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
    

class PerfilPotencia:
    segundos_por_mes=config.SEGUNDOS_POR_MES
    def __init__(self, perfil_path=config.perfil_path):
        self.perfil_path = perfil_path

    def leitura_csv(self):
        dados = pd.read_csv(self.perfil_path)
        Pot_out_perfil = dados['TractionPower'].values
        Pot_in_perfil = dados['BrakingPower'].values
        Pot_total = Pot_out_perfil - Pot_in_perfil #forma pyBaMM
        
        return 1000*Pot_total

    def carregar_perfil_potencia(self, segundos_totais=None):
        if segundos_totais == None:
            segundos_totais = PerfilPotencia.segundos_por_mes
        
        Pot_total = -self.leitura_csv()
        tam_raw = len(Pot_total)

        n_repeticoes = int(np.ceil(segundos_totais / tam_raw))
        #1683 vezes que o perfil de potência é replicado 
        Pot_total_expandido = np.tile(Pot_total, n_repeticoes)[:segundos_totais]
        mission_day_pot = Pot_total_expandido 
        
        return mission_day_pot

    def calcular_energia_tracao_mensal(self):
        """
        Energia de tracao total exigida pelas rodas no mes (MWh), independente
        de Np - nao passa pelo BESS/PLECS. Usada como base do custo Cd(t)
        (diesel sem hibridizacao). O perfil se repete identico todo mes, entao
        o valor e o mesmo em todos os meses do horizonte.
        """
        dados = pd.read_csv(self.perfil_path)
        Pot_out_perfil = dados['TractionPower'].values
        tam_raw = len(Pot_out_perfil)

        n_repeticoes = int(np.ceil(PerfilPotencia.segundos_por_mes / tam_raw))
        Pot_out_expandido = np.tile(Pot_out_perfil, n_repeticoes)[:PerfilPotencia.segundos_por_mes]

        return np.sum(Pot_out_expandido) / 3600 / 1000  # MWh


'''try:
    from rich import inspect
    inspect(Bateria)
except: pass'''






