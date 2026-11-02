import config
from SOCestimation_testePOO import Simulacao, PerfilPotencia
try:
    from rich import print
    from rich.table import Table
    RICH_DISPONIVEL = True
except ImportError:
    RICH_DISPONIVEL = False


class LCC:
    """
    Calculo do custo de ciclo de vida (LCC) do BESS em funcao de Np, a partir
    dos resultados mensais de SOCestimation_teste.simular_missao.
    """
    def __init__(
        self,
        n_meses=config.n_meses,
        Np = config.Npb_padrao,
        perfil_path=config.perfil_path,
        verbose = False
    ):
        
        self.n_meses = n_meses
        Perfil_pot = PerfilPotencia(perfil_path=perfil_path)
        self.E_tracao_MWh = Perfil_pot.calcular_energia_tracao_mensal()
        mission_day_pot = Perfil_pot.carregar_perfil_potencia()

        self.Np = round(Np)
        self.im = (1 + config.Y) ** (1 / 12) - 1

        self.simulado = Simulacao(Npb=self.Np, mission_day_pot=mission_day_pot, 
                                  n_meses=self.n_meses)
        self.simulado.executar(verbose=verbose)


    def capex_pack(self):
        """
        Retorna (CAPEX total, custo do pack isolado) em USD, para um dado Np.
        CAPEX = pack ($/kWh) + conversor ($/kW); o custo do pack isolado e
        reaproveitado tambem para troca endogena e valor residual.
        """
        einstalada_kWh = config.Nsb * self.Np * config.Vcell * config.CAh0 / 1000
        pmax_kW = einstalada_kWh * config.Cr
        custo_pack = config.Cpack * einstalada_kWh
        custo_conversor = config.Cconv * pmax_kW

        return custo_pack + custo_conversor, custo_pack


    def iniciando_variaveis_custo(self):
        self.CAPEX_Np, self.custo_pack = self.capex_pack()
        self.Coper_mensal = config.Cd_price * self.E_tracao_MWh * 1000  # MWh -> kWh
        self.Coper_total_vp = 0.0
        self.Cer_total_vp = 0.0
        self.custo_trocas_vp = 0.0
        self.custo_total = self.CAPEX_Np
        
    
    def mostra_valores(self, Cresid_vp, Eregen_total, Ethroughput_total, custo_diesel_liquido_vp, Cf_T):
        if RICH_DISPONIVEL:
            tabela = Table(title=f"LCC — Np = {self.Np:g}")

            tabela.add_column("Indicador", justify="left")
            tabela.add_column("Valor", justify="right")
            tabela.add_column("Unidade", justify="left")

            tabela.add_row("Custo total f(Np)",f"{self.custo_total:,.2f}","USD")
            tabela.add_row("CAPEX",f"{self.CAPEX_Np:,.2f}","USD")
            tabela.add_row("Valor residual (VP)",f"{Cresid_vp:,.2f}","USD")
            tabela.add_row("Custo diesel líquido (VP)",f"{custo_diesel_liquido_vp:,.2f}","USD")
            tabela.add_row("Custo das trocas (VP)",f"{self.custo_trocas_vp:,.2f}","USD")
            tabela.add_row("Receita de regeneração (VP)",f"{self.Cer_total_vp:,.2f}","USD")
            tabela.add_row("Energia regenerada",f"{Eregen_total:.6g}","MWh")
            tabela.add_row("Energy throughput",f"{Ethroughput_total:.6g}","MWh")
            tabela.add_row("Capacity fade final", f"{Cf_T:.4g}","%")
            print(tabela)
            print(f"Horizonte: {self.n_meses}")

        else: 
            print(f"[LCC] Np={self.Np:g}")
            print(f"  f(Np)={self.custo_total:,.2f} USD | CAPEX={self.CAPEX_Np:,.2f} USD | "
            f"Cresid(VP)={Cresid_vp:,.2f} USD | Cf_final={Cf_T:.4g}%")
            print(f"  Eregen_total={Eregen_total:.6g} MWh | Custo_Eregen(CER, VP)={self.Cer_total_vp:,.2f} USD | "
                    f"Ethroughput_total={Ethroughput_total:.6g} MWh")
            print(f"  Custo_diesel_liquido(VP)={custo_diesel_liquido_vp:,.2f} USD "
                    f"(bruto VP={self.Coper_total_vp:,.2f} - CER VP={self.Cer_total_vp:,.2f})")
            print(f"custo_trocas_total(VP)={self.custo_trocas_vp:,.2f} USD")
            
        
    
    def custo_ciclo_vida(self):
        """
        f(Np): custo de ciclo de vida do BESS, somando CAPEX, fluxo mensal
        descontado (diesel liquido + trocas endogenas) e subtraindo o valor
        residual do pack ao final do horizonte.

        Np e o numero de strings em paralelo - fisicamente precisa ser inteiro.
        Arredondado aqui (nao no chamador) porque o Brent, sendo um otimizador
        continuo, testa valores fracionarios; este e o unico ponto por onde
        todo Np passa antes de tocar a simulacao/CAPEX.

        perfil_cache (PerfilCache, opcional) evita reler o Excel do perfil e
        reescrever perfil_potencia.csv a cada avaliacao de Np - gerado uma vez
        por preparar_perfil(perfil_path) e reaproveitado pelo grid search e
        pelo Brent num mesmo otimizar(). Validado contra perfil_path recebido
        aqui: se o cache foi gerado para outro perfil, levanta erro em vez de
        silenciosamente simular a mina errada (perfil_potencia.csv em disco
        nao seria o mesmo do perfil_path pedido).
        """
        self.iniciando_variaveis_custo()
        
        for t in range(1, self.n_meses + 1):
            i = t - 1
            fator_desc = 1 / (1 + self.im) ** t

            Cer_t = config.CER_price * self.simulado.Eregen_aproveitada[i] * 1000  # MWh -> kWh
            troca_t = self.custo_pack if self.simulado.troca[i] else 0.0

            self.Coper_total_vp += self.Coper_mensal * fator_desc
            
            self.Cer_total_vp += Cer_t * fator_desc
            self.custo_trocas_vp += troca_t * fator_desc

            fluxo_mes = self.Coper_mensal - Cer_t + troca_t
            self.custo_total += fluxo_mes * fator_desc

        Cf_T = self.simulado.Cf[-1]
        Ethroughput_total, Eregen_total = self.simulado.calcular_energia_total()

        Cresid = self.custo_pack * (1 - Cf_T / config.Cf_EOL)
        Cresid_vp = Cresid / (1 + self.im) ** self.n_meses
        self.custo_total -= Cresid_vp
        custo_diesel_liquido_vp = self.Coper_total_vp - self.Cer_total_vp

        self.mostra_valores(Cresid_vp, Eregen_total, Ethroughput_total, custo_diesel_liquido_vp, Cf_T)
        
        return self.custo_total


if __name__ == "__main__":
    lcc = LCC(Np=config.Npb_padrao)
    custo = lcc.custo_ciclo_vida()
    print(f"\nf(Np={config.Npb_padrao:g}) = {custo:,.2f} USD")
