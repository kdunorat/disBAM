import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from math import floor


class DepthAnalyzer:
    """Análises da profundidade do sequenciamento"""
    def __init__(self, df_depth, filename):
        self.df_depth = df_depth
        self.filename = filename

    def _th_cross(self):
        """Vezes em que houve uma queda da região acima dos 400 para abaixo dos 400"""
        low = self.df_depth['depth'].gt(400)
        crosses = low[low.ne(low.shift(-1))][:-1]
        return f'crossed the threshold line of depth 400: {floor(crosses.count() / 2)} times'

    def _reads_below(self):
        """Porcentagem de reads abaixo dos 400 de profundidade"""
        indexes = np.where(self.df_depth['depth'] < 400)
        percentage = (len(indexes[0]) * 100) / len(self.df_depth)

        return int(percentage)

    def _generate_graph(self):
        threshold = 7000

        chart = sns.relplot(data=self.df_depth, kind='line', height=10, aspect=2,
                            y='depth', x='position')

        chart.axes[0][0].axhline(y=self.df_depth.depth.mean(), color='red', linewidth=2, ls=':')
        chart.axes[0][0].axhline(y=threshold, color='yellow', linewidth=2, alpha=.7)

        plt.savefig(f'output/{self.filename}_depth.png')

    def get_final_text(self):
        cross = self._th_cross()
        percentage = self._reads_below()
        return f'{cross}\t{percentage}%\n'

    def run_depth_analysis(self):
        self._generate_graph()

        return self.get_final_text()
