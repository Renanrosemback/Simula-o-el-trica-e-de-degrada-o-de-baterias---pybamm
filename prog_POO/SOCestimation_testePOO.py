import time
import numpy as np

import config
from battery_circuit_POO import Bateria
#from bibliotecapyBaMM import Bateria #mudei a classe bateria
from battery_circuit_POO import PerfilPotencia
from lifetime_degradationPOO import BatteryDegradation


try:
    from rich import print
    from rich.panel import Panel
    from rich.columns import Columns
    from rich.console import Group
    RICH_DISPONIVEL = True
except ImportError:
    RICH_DISPONIVEL = False





class Simulacao:
    def __init__(
        self,
        Npb=config.Npb_padrao,
        n_meses=config.n_meses,
        perfil_path=config.perfil_path,
        Tbat=config.Tbat,
        Cf_replace=config.Cf_EOL,
        mission_day_pot=None
    ):

        self.Npb = Npb
        self.n_meses = n_meses
        self.Perfil_pot = PerfilPotencia(perfil_path)
        self.Tbat = Tbat
        self.Cf_replace = Cf_replace

        # Perfil de potência
        if mission_day_pot is None:
            self.mission_day_pot = self.Perfil_pot.carregar_perfil_potencia()
        else:
            self.mission_day_pot = mission_day_pot
            

        # Objetos principais da simulação
        self.bateria = Bateria(Npb=Npb)

        self.degradacao = BatteryDegradation(
            Tbat=Tbat,
            C_fade_cyc_in=0.0,
            C_fade_Cal_in=0.0
        )

        # Histórico
        self._inicializar_historico()
        
    def _mostrar_estado(self, i=0, inicial=False, final=False):
        if inicial:
            print('[green]Começando simulação com pack específicado:[/]')
            conteudo = f" Npb={self.Npb:g} "
            conteudo += f'Rs={self.bateria.Rs:.6g} '
            conteudo += f"R1={self.bateria.R1:.6g} Q1={self.bateria.Q1:.6g}"
            conteudo +=  f"R2={self.bateria.R2:.6g} Q2={self.bateria.Q2:.6g}"
            conteudo +=  f"CAh={self.bateria.CAh:.6g} SOCi={self.bateria.soc:.4f} "
            if RICH_DISPONIVEL:
                painel = Panel(conteudo, title= 'Dados iniciais')
                print(painel)
            else:
                print('Baixe a biblioteca rich para melhorar a qualida do print')
                print(conteudo)

        if final:
            print(f"Tempo total simulacao={self.tempo_total_simulacao:.2f}s")
                        
            if len(self.meses_troca) > 0:
                print(f"\nTrocas de pack nos meses: {self.meses_troca}")
            else:
                print("\nNenhuma troca de pack ao longo do horizonte simulado.")
            print(f"Total: {len(self.meses_troca)} trocas. ") 
            print("[bold green]Processo finalizado.[/]")
        
        elif final == False and inicial == False:
            conteudo_desgaste = f"troca={self.troca[i]} \n"
            conteudo_desgaste += f"\n -> Cf={self.Cf[i]:.6g}% \n"
            conteudo_desgaste += f"C_cal={self.C_cal[i]:.6g}% \n"
            conteudo_desgaste += f"C_cyc={self.C_cyc[i]:.6g} \n"
            
            conteudo_estado = f'Rs={self.bateria.Rs:.6g} \n'
            conteudo_estado += f"R1={self.bateria.R1:.6g} Q1={self.bateria.Q1:.6g} \n"
            conteudo_estado +=  f"R2={self.bateria.R2:.6g} Q2={self.bateria.Q2:.6g} \n"
            conteudo_estado +=  f"CAh={self.bateria.CAh:.6g}\n "
            
            conteudo_potencia = f"Ethroughput={self.Ethroughput[i]:.6g} MWh \n"
            conteudo_potencia+= f"Eregen={self.Eregen_aproveitada[i]:.6g} MWh \n"
            conteudo_potencia += f"SOC_final={self.SOC_final[i]:.6g} \n"
        
            conteudo_tempo = f"Tempo gasto execução mês {i+1} = {self.tempo_mes[i]:.2f}s \n"
            
            if RICH_DISPONIVEL:
                painel_desgaste = Panel(conteudo_desgaste, title='[bold red]Desgaste[/]')
                painel_estado = Panel(conteudo_estado, title='[bold blue]Estado bateria[/]')
                painel_potencia = Panel(conteudo_potencia, title='[bold green]Potência[/]')
                colunas = Columns([painel_desgaste, painel_estado, painel_potencia])
                grupo = Group(conteudo_tempo, colunas)
                
                painel_principal = Panel(
                        grupo, 
                        title=f'Mes {i+1}/{self.n_meses}', 
                        expand=False
                    )
                print(painel_principal)
            else:
                print(conteudo_desgaste, conteudo_estado, conteudo_potencia)
                      
    def _inicializar_historico(self):
        n = self.n_meses
        self.Cf = np.zeros(n)
        self.C_cal = np.zeros(n)
        self.C_cyc = np.zeros(n)
        
        self.Ethroughput = np.zeros(n)
        self.Eregen_aproveitada = np.zeros(n)
        
        self.troca = np.zeros(n, dtype=bool)
        self.meses_troca = []
        
        self.SOC_final = np.zeros(n)
        
        self.Rs_hist = np.zeros(n)
        self.R1_hist = np.zeros(n)
        self.R2_hist = np.zeros(n)
        self.Q1_hist = np.zeros(n)
        self.Q2_hist = np.zeros(n)
        self.CAh_hist = np.zeros(n)
        
        self.tempo_engine = np.zeros(n)
        self.tempo_degradacao = np.zeros(n)
        self.tempo_resultados = np.zeros(n)
        self.tempo_mes = np.zeros(n) #so o que estou utilizando por enquanto é esse
        
    def _calcular_resultados_energeticos(self, i, P_bess):
        self.Ethroughput[i] = \
            np.sum(np.abs(P_bess)) / 3600 / 1e6
    
        self.Eregen_aproveitada[i] = \
            np.sum(P_bess[P_bess > 0]) / 3600 / 1e6
            
    def _processar_degradacao(self, i, resultado_degradacao):
        self.C_cal[i] = resultado_degradacao['C_fade_cal']
        self.C_cyc[i] = resultado_degradacao['C_fade_cyc']
        self.Cf[i] = self.C_cyc[i] + self.C_cal[i]

    def _deve_trocar(self, i):
        return self.Cf[i] >= self.Cf_replace
            
       
    def _realizar_troca(self, i):
        self.troca[i] = True
        self.meses_troca.append(i)
        self.bateria.resetar()
        self.degradacao.resetar()
        
    def _atualizar_bateria(self, resultado_degradacao):
        self.bateria.aplicar_degradacao(**resultado_degradacao)
        self.degradacao._atualizar_acumulado(
            resultado_degradacao['C_fade_cyc'],
            resultado_degradacao['C_fade_cal'])
        
    def _guardar_estado(self, i):
        self.Rs_hist[i] = self.bateria.Rs
        self.R1_hist[i] = self.bateria.R1 
        self.R2_hist[i] = self.bateria.R2
        self.Q1_hist[i] = self.bateria.Q1
        self.Q2_hist[i] = self.bateria.Q2
        self.CAh_hist[i] = self.bateria.CAh

    def calcular_energia_total(self):
        Eregen_total = float(np.sum(self.Eregen_aproveitada))
        Ethroughput_total = float(np.sum(self.Ethroughput))  
        return Ethroughput_total, Eregen_total
    
         
    def executar(self, verbose=True):
        inicio = time.perf_counter()
        if verbose:
            self._mostrar_estado(inicial=True) 
    
        for i in range(self.n_meses):
    
            inicio_mes = time.perf_counter()
    
            # 1 - FAZENDO SIMUAÇÃO BATERIA:
            SOC_pct, P_bess, _, _ = self.bateria.simular_celula(self.mission_day_pot)
        
            # 2 - RESULTADOS ENERGÉTICOS
            mission_day_SOC = np.asarray(SOC_pct)
            P_bess = np.asarray(P_bess)
            self._calcular_resultados_energeticos(i, P_bess)
            
            # 3 - DEGRADAÇÃO
            resultado_degradacao = self.degradacao.lifetime_degradation(mission_day_SOC)
    
            # 4 - PROCESSAR DEGRADAÇÃO
            self._processar_degradacao(i, resultado_degradacao)
            self.SOC_final[i] = mission_day_SOC[-1]
    
            # 5 - VERIFICAR TROCA
            if self._deve_trocar(i): 
                self._realizar_troca(i)
    
            else: 
                self._atualizar_bateria(resultado_degradacao)
    
            # 6 - GUARDAR ESTADO
            self._guardar_estado(i)
    
            self.tempo_mes[i] = \
                time.perf_counter() - inicio_mes
                
            if verbose:
                self._mostrar_estado(i) 
    
        self.tempo_total_simulacao = \
            time.perf_counter() - inicio
        
        if verbose:
            self._mostrar_estado(final=True)

            
    
        


if __name__ == "__main__":
    resultado = Simulacao(n_meses=10)
    resultado.executar()

    



