"""
Constantes economicas e fisicas para o calculo de custo de ciclo de vida
(LCC) do dimensionamento do BESS. Fontes e decisoes documentadas na
memoria do projeto (project_bess_lcc_optimization).
"""

# Custos de investimento
Cpack = 265.0   # $/kWh (BNEF/Energy Innovation 2020, HDV, ja inclui BOS)
Cconv = 890.0   # $/kW (conversor bidirecional, proxy estacionario, arXiv 2107.01032)

# Custo de energia (diesel-equivalente entregue as rodas)
Cd_price = 0.70   # $/kWh
CER_price = 0.70  # $/kWh (mesmo valor de Cd_price)

# Degradacao / fim de vida
Cf_EOL = 20.0    # % de fade de capacidade para troca/fim de vida

# Taxa de desconto - PLACEHOLDER, confirmar taxa real antes da rodada final
Y = 0.1425         

# Limites de busca de Np - PLACEHOLDER exploratorio, sem restricao fisica definida ainda
Np_bounds = (1, 41000)  
n_grid = 6

# Horizonte de simulacao usado no calculo do LCC (meses). 4 por enquanto,
# so pra validar a integracao SOCestimation -> lcc -> otimizacao antes de
# ir para o horizonte final de 60 meses.
n_meses = 36 #lembrar de voltar para 36

# Perfil de potencia da mina usado na simulacao (uma mina por rodada - ver
# memoria do projeto). Trocar aqui para rodar outra mina - CSVs padronizados
# (colunas TractionPower/BrakingPower) disponiveis nesta pasta:
#   perfil_1.csv = AGREGADOR ANALYSIS
#   perfil_2.csv = CR-3112 (mina usada em toda a validacao contra a PLECS)
#   perfil_3.csv = UMAX
perfil_path = 'perfil_2.csv'

# Duracao do mes usada tanto para expandir o perfil de potencia
# (SOCestimation_teste.carregar_perfil_potencia) quanto para calcular a
# energia de tracao mensal (lcc.calcular_energia_tracao_mensal) - fonte
# unica para as duas funcoes nao divergirem silenciosamente se um dia for
# alterada em um lugar so.
SEGUNDOS_POR_MES = 30 * 24 * 60 * 60  # 2_592_000

# Np usado nos blocos "if __name__ == '__main__'" de smoke-test
# (SOCestimation.py/SOCestimation_teste.py/lcc.py rodados isolados) - a
# otimizacao de verdade busca o Np sozinha via otimizacao.otimizar().
Npb_padrao = 2600.0

# Parametros fisicos do pack (mesmos valores usados como default em
# SOCestimation_teste.simular_missao - passados explicitamente ao chamar,
# pra nao depender de defaults ficarem sincronizados por coincidencia)
Nsb = 390.0     # Ns fixo em 1
Vcell = 3.3   # V, tensao nominal da celula
CAh0 = 2.5    # Ah, capacidade nominal por celula
Cr = 4.0      # C-rate maximo da celula
Tbat = 35.0   # C, temperatura da bateria (fixa, sem modelo termico ainda)

# Janela de operacao do SOC e SOC inicial da missao (%)
SOCmin = 20.0
SOCmax = 80.0
SOCi0 = 20.0

# Parametros eletricos de fabrica do pack (Rs em Ohm, R1/R2 em Ohm, Q1/Q2 em
# Ah - entram em calc_zarc_params para montar os ramos ZARC1/ZARC2). Sao os
# valores de "pack novo", usados tanto na primeira simulacao quanto para
# resetar o estado a cada troca endogena.
Rs0 = 0.017
R10 = 0.003
R20 = 0.06
Q10 = 2.8
Q20 = 890.0

# Curva de tensao de circuito aberto OCV(SOC) da celula - fixa, caracteristica
# da quimica/celula usada (nao depende de Np/mina). SOC em fracao (0 a 1),
# igualmente espacada por SOC_STEP - permite interpolacao O(1) por indice em
# battery_circuit._ocv em vez de busca a cada um dos ~2.6 milhoes de passos
# por mes.
OCV_TABLE = (2.5, 3.25, 3.45, 3.49, 3.5, 3.53, 3.54, 3.55, 3.56, 3.57, 3.571,
             3.573, 3.574, 3.576, 3.577, 3.578, 3.58, 3.58, 3.582, 3.585, 3.62)
SOC_STEP = 0.05

# Passo de discretizacao (s) do motor eletrico (battery_circuit.py)
DT = 1.0
