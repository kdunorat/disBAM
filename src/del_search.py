import re


class DelSearch:
    """Procura prováveis deleções a partir de regiões vizinhas que sofreram softclip"""
    def __init__(self, df, data):
        self.df = df
        self.data = data
        self.depth_set = set()
        self.analyzed_regions = set()
        self.text = 'Region\n'

    def _get_depth(self):
        """Obtém apenas as regiões de softclip com profundidade maior ou igual a 10"""
        depth_range = self.df['fmap'].value_counts()
        depth_range = depth_range.where(depth_range >= 20)
        depth_range = depth_range.dropna()
        for key in depth_range.keys():
            self.depth_set.add(key)

    @staticmethod
    def _get_soft_seq(cigar, seq):
        """Gera 'strings' com a região inicial e final de softclip e também com a última região mapeada"""
        pattern = r'(\d+)(S|M)'
        soft_qt = re.findall(pattern, cigar)
        first_soft, last_soft, mapped = False, False, False
        for tupla in soft_qt:
            if tupla[1] == 'M':
                mapped = seq[:int(tupla[0])]
                seq = seq[int(tupla[0]):]
            elif tupla[1] == 'S' and mapped:
                if int(tupla[0]) > last_soft:
                    last_soft = seq[:int(tupla[0])]
            else:
                first_soft = seq[:int(tupla[0])]
                seq = seq[int(tupla[0]):]

        return first_soft, last_soft, mapped

    def _analyze_reads(self, row, last_soft, mapped, count):
        """Analiza par a par em busca de provaveis deleções"""
        distance = 0
        for row_2 in self.data[count:]:
            if row == row_2:
                continue
            first_soft_2, _, mapped_2 = self._get_soft_seq(row_2['cigar'], row_2['seq'])
            if first_soft_2:
                initpos_1 = row['fmap']
                initpos_2 = row_2['pos']
                if initpos_2 > initpos_1:
                    distance = initpos_2 - initpos_1
                if distance != 0 and distance < 800:
                    if last_soft in mapped_2 and first_soft_2 in mapped:
                        if initpos_1 not in self.analyzed_regions:
                            prompt = f'{initpos_1}>{initpos_2}\n'
                            self.text += prompt
                            self.analyzed_regions.add(initpos_1)

    def _create_log(self):
        """Cria o arquivo final"""
        with open('output/log.txt', 'w') as f:
            f.write(self.text)

    def run_analysis(self):
        """Roda o processo inteiro"""
        count = 0
        self._get_depth()
        for row in self.data:
            if row['fmap'] not in self.depth_set:
                count += 1
                continue
            _, last_soft, mapped = self._get_soft_seq(row['cigar'], row['seq'])
            if last_soft:
                self._analyze_reads(row, last_soft, mapped, count)
            count += 1
        self._create_log()
