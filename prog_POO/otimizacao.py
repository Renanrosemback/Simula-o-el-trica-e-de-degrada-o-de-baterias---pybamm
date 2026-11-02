from datetime import datetime

import numpy as np
from scipy.optimize import minimize_scalar

import config
from lcc import LCC


class Otimizacao:

    def __init__(
        self,
        n_meses=config.n_meses,
        n_grid=config.n_grid,
        perfil_path=config.perfil_path,
        Np_bounds=config.Np_bounds,
        verbose=True
    ):

        self.n_meses = n_meses
        self.n_grid = n_grid
        self.perfil_path = perfil_path
        self.Np_min, self.Np_max = Np_bounds
        self.verbose = verbose

        self.grid = None
        self.custos_grid = None

        self.Np_melhor_grid = None
        self.custo_melhor_grid = None
        self.bounds_refino = None

        self.resultado_brent = None
        self.Np_otimo = None
        self.custo_otimo = None


    def _imprimir_cabecalho(self):

        print("=" * 70)
        print(
            f"[otimizacao] Inicio: "
            f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
        )
        print(f"[otimizacao] Perfil: {self.perfil_path}")
        print(f"[otimizacao] Tbat: {config.Tbat:g} °C")
        print(f"[otimizacao] Np: {self.Np_min:g} → {self.Np_max:g}")
        print(f"[otimizacao] Grid: {self.n_grid} pontos")
        print(f"[otimizacao] Horizonte: {self.n_meses} meses")
        print("=" * 70)


    def _avaliar_Np(self, Np):

        lcc = LCC(
            n_meses=self.n_meses,
            Np=Np,
            perfil_path=self.perfil_path
        )

        return lcc.custo_ciclo_vida()


    def _criar_grid(self):

        self.grid = np.unique(
            np.round(
                np.linspace(
                    self.Np_min,
                    self.Np_max,
                    self.n_grid
                )
            ).astype(int)
        )


    def _avaliar_grid(self):

        print(
            f"[otimizacao] Grid inicial "
            f"({len(self.grid)} pontos):"
        )

        self.custos_grid = np.array([
            self._avaliar_Np(Np)
            for Np in self.grid
        ])


    def _encontrar_melhor_grid(self):

        i_melhor = int(np.argmin(self.custos_grid))

        self.Np_melhor_grid = self.grid[i_melhor]
        self.custo_melhor_grid = self.custos_grid[i_melhor]

        passo = self.grid[1] - self.grid[0]

        self.bounds_refino = (
            max(
                self.Np_min,
                self.Np_melhor_grid - passo
            ),
            min(
                self.Np_max,
                self.Np_melhor_grid + passo
            )
        )

        print(
            f"\n[otimizacao] Melhor ponto do grid: "
            f"Np={self.Np_melhor_grid} "
            f"f={self.custo_melhor_grid:,.2f} USD"
        )


    def _refinar_brent(self):

        print(
            f"[otimizacao] Refinando com Brent em "
            f"[{self.bounds_refino[0]:g}, "
            f"{self.bounds_refino[1]:g}]..."
        )

        self.resultado_brent = minimize_scalar(
            self._avaliar_Np,
            bounds=self.bounds_refino,
            method="bounded",
            options={"xatol": 10}
        )

        self.Np_otimo = round(self.resultado_brent.x)
        self.custo_otimo = self.resultado_brent.fun


    def resultado(self):

        return {
            "grid": self.grid,
            "custos_grid": self.custos_grid,
            "Np_melhor_grid": self.Np_melhor_grid,
            "custo_melhor_grid": self.custo_melhor_grid,
            "bounds_refino": self.bounds_refino,
            "Np_otimo": self.Np_otimo,
            "custo_otimo": self.custo_otimo,
            "resultado_brent": self.resultado_brent,
        }


    def otimizar(self):

        self._imprimir_cabecalho()

        self._criar_grid()
        self._avaliar_grid()
        self._encontrar_melhor_grid()
        self._refinar_brent()

        print(
            f"\n[otimizacao] Np otimo: "
            f"{self.Np_otimo} "
            f"f(Np) = {self.custo_otimo:,.2f} USD"
        )

        return self.resultado()


if __name__ == "__main__":

    otimizador = Otimizacao()
    resultado = otimizador.otimizar()