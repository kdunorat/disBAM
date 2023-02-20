import re


class DelSearch:
    """Procura prováveis deleções a partir de regiões vizinhas que sofreram softclip"""
    def __init__(self, df, data):
        self.df = df
        self.data = data
        self.depth_mapped = set()
        self.analyzed_regions = set()
        self.text = 'Region\tSoft Cliped\n'

    def _get_depth_mapped(self):
        """Obtém apenas as regiões de softclip com profundidade maior ou igual a 10"""
        depth_range = self.df['fmap'].value_counts()
        depth_range = depth_range.where(depth_range >= 20)
        depth_range = depth_range.dropna()
        for key in depth_range.keys():
            self.depth_mapped.add(key)

    def _get_depth_soft(self, fmap):
        new = self.df[self.df['fmap'].isin([fmap])]
        depth_range = new['fpos'].value_counts()
        depth_soft = depth_range.keys()
        depth_soft = depth_soft[0]
        return depth_soft

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

    def _analyze_reads(self, read, last_soft, mapped, count):
        """Analiza par a par em busca de provaveis deleções"""
        distance = 0
        for read_2 in self.data[count:]:
            first_soft_2, _, mapped_2 = self._get_soft_seq(read_2['cigar'], read_2['seq'])
            if first_soft_2:
                initpos_1 = read['fmap']
                initpos_2 = read_2['pos']
                if initpos_2 > initpos_1:
                    distance = initpos_2 - initpos_1
                if distance != 0 and distance < 180:
                    if last_soft in mapped_2 and first_soft_2 in mapped:
                        if initpos_1 not in self.analyzed_regions:
                            prompt = f'{initpos_1}>{initpos_2}\t{last_soft}\n'
                            self.text += prompt
                            self.analyzed_regions.add(initpos_1)

    def _create_log(self):
        """Cria o arquivo final"""
        with open('output/log.txt', 'w') as f:
            f.write(self.text)

    def run_analysis(self):
        """Roda o processo inteiro"""
        count = 1
        self._get_depth_mapped()
        for read in self.data:
            if read['fmap'] not in self.depth_mapped:
                count += 1
                continue
            _, last_soft, mapped = self._get_soft_seq(read['cigar'], read['seq'])
            if last_soft:
                depth_soft = self._get_depth_soft(read['fmap'])
                if read['fpos'] == depth_soft:
                    self._analyze_reads(read, last_soft, mapped, count)
            count += 1
        self._create_log()
