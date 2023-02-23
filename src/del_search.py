from process_data import get_soft_seq


class DelSearch:
    """Procura prováveis deleções a partir de regiões vizinhas que sofreram softclip"""
    def __init__(self, df, data):
        self.df = df
        self.data = data
        self.depth_set_mapped, self.depth_set_pos = set(), set()
        self.analyzed_regions = set()
        self.text = 'Region\tSoft Clip\n'

    def _get_depth_mapped(self):
        """Obtém apenas as regiões de softclip com profundidade maior ou igual a 10"""
        depth_mapped = self.df['fmap'].value_counts()
        depth_pos = self.df['pos'].value_counts()
        depth_mapped = depth_mapped.where(depth_mapped >= 20)
        depth_pos = depth_pos.where(depth_pos >= 20)
        depth_mapped = depth_mapped.dropna()
        depth_pos = depth_pos.dropna()
        for key in depth_mapped.keys():
            self.depth_set_mapped.add(key)
        for key in depth_pos.keys():
            self.depth_set_pos.add(key)

    def _get_depth_soft(self, fmap):
        new = self.df[self.df['fmap'].isin([fmap])]
        depth_soft = new['fpos'].value_counts()
        depth_soft = depth_soft.keys()
        depth_soft = depth_soft[0]
        return depth_soft

    def _analyze_reads(self, read, last_soft, mapped, count):
        """Analiza par a par em busca de provaveis deleções"""
        distance = 0
        for read_2 in self.data[count:]:
            if read_2['pos'] not in self.depth_set_pos:
                continue

            first_soft_2, _, mapped_2 = get_soft_seq(read_2['cigar'], read_2['seq'])
            if first_soft_2 and len(first_soft_2) > 50:
                fmap_left = read['fmap']
                pos_right = read_2['pos']
                if pos_right > fmap_left:
                    distance = pos_right - fmap_left
                if distance != 0 and distance < 350:
                    self._match_check(mapped, first_soft_2, mapped_2,
                                      fmap_left, pos_right, last_soft)

    def _match_check(self, mapped, first_soft_2, mapped_2, fmap_left, pos_right, last_soft):
        if last_soft in mapped_2 or first_soft_2 in mapped:
            if fmap_left not in self.analyzed_regions:
                prompt = f'{fmap_left}>{pos_right}\t{last_soft}\n'
                self.text += prompt
                self.analyzed_regions.add(fmap_left)

    def _create_log(self):
        """Cria o arquivo final"""
        if len(self.text) <= 17:
            self.text += '----Não foram encontradas deleções extras----'
        with open('src/output/log.txt', 'w') as f:
            f.write(self.text)

    def run_analysis(self):
        """Roda o processo inteiro"""
        count = 1
        self._get_depth_mapped()
        for read in self.data:
            if read['fmap'] not in self.depth_set_mapped:
                count += 1
                continue
            _, last_soft, mapped = get_soft_seq(read['cigar'], read['seq'])
            if last_soft and len(last_soft) > 50:
                depth_soft = self._get_depth_soft(read['fmap'])
                if read['fpos'] == depth_soft:
                    self._analyze_reads(read, last_soft, mapped, count)
            count += 1
        self._create_log()
