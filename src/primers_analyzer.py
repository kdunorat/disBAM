import subprocess
import pandas as pd


class PrimersAnalyzer:
    def __init__(self, file, absolute_path):
        self.absolute_path = absolute_path
        self.file = file

    def _create_stats(self):
        subprocess.run(['bash', f'{self.absolute_path}/src/sh_scripts/primerstats.sh', self.file])

    def _extract_from_text(self):
        self.amp, self.frperc, self.fdepth, self.freads, self.fpcov, self.famp = [], [], [], [], [], []
        with open(f'{self.absolute_path}/src/primers/stats.txt', 'r') as stats:
            for line in stats:
                if line.startswith("AMPLICON"):
                    self.amp.append(line.strip())
                elif line.startswith("FRPERC"):
                    self.frperc = self._process_line(line)
                elif line.startswith("FDEPTH"):
                    self.fdepth = self._process_line(line)
                elif line.startswith("FREADS"):
                    self.freads = self._process_line(line)
                elif line.startswith("FPCOV"):
                    self.fpcov = self._process_line(line)
                elif line.startswith("FAMP"):
                    self.famp.append(self._process_line(line))
            self.famp.pop(0)

    @staticmethod
    def _process_line(line):
        line_list = line.split()
        line_list = line_list[2:]
        return line_list
    
    def _create_data(self):
        primers_features = pd.DataFrame([linha.split("\t") for linha in self.amp])
        primers_features = primers_features.iloc[:, 3:]
        primers_features.rename(columns={3: 'left', 4: 'right'}, inplace=True)

        # Inserting other features
        primers_features['frperc'] = self.frperc
        primers_features['fdepth'] = self.fdepth
        primers_features['freads'] = self.freads
        primers_features['fpcov'] = self.fpcov
        primers_features['famp'] = self.famp

    def run(self):
        # self._create_stats()
        self._extract_from_text()


if __name__ == '__main__':
    prim = PrimersAnalyzer('aln.virus.sorted_NC_045512.2.bam', '/home/kdunorat/PycharmProjects/disBAM')
    prim.run()



