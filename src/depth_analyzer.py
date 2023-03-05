import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from math import floor


class DepthAnalyzer:
    def __init__(self, df_depth):
        self.df_depth = df_depth

    def _th_cross(self, threshold):
        low = self.df_depth['depth'].gt(threshold)
        crosses = low[low.ne(low.shift(-1))][:-1]
        print(f'crossed the threshold line of {threshold} reads {floor(crosses.count() / 2)} times')

    def _reads_below(self, threshold):
        """Porcentagem de reads abaixo dos 400 de profundidade"""
        indexes = np.where(self.df_depth['depth'] < threshold)
        percentage = (len(indexes[0]) * 100) / len(self.df_depth)

        return percentage

    def generate_graph(self, filename):
        threshold = 7000

        chart = sns.relplot(data=self.df_depth, kind='line', height=10, aspect=2,
                            y='depth', x='position')

        chart.axes[0][0].axhline(y=self.df_depth.depth.mean(), color='red', linewidth=2, ls=':')
        chart.axes[0][0].axhline(y=threshold, color='yellow', linewidth=2, alpha=.7)

        plt.savefig(f'output/{filename}_depth.png', transparent=True)
