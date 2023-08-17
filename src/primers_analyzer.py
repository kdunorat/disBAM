import subprocess
import pandas as pd


class PrimersAnalyzer:
    def __init__(self, file, absolute_path):
        self.absolute_path = absolute_path
        self.file = file

    def _create_stats(self):
        subprocess.call(['bash', f'{self.absolute_path}/sh_scripts/primers_analyzer.sh', self.file])

    @staticmethod
    def _create_dataframe():
        with open('primers/stats.txt', 'r') as f:
            stats = f.read()
        print(stats)

    def run(self):
        self._create_stats()
        #self._create_dataframe()


if __name__ == '__main__':
    prim = PrimersAnalyzer('aln.virus.sorted_NC_045512.2.bam', '/home/kdunorat/PycharmProjects/disBAM')
    prim.run()



