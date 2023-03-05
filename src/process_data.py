import subprocess
import pandas as pd
import re
import os


class InputHandler:
    """Processa os inputs necessários ao funcionamento do programa"""
    def __init__(self, file, absolute_path):
        self.filename = file
        self.absolute_path = absolute_path

    def create_data(self):
        """Cria um dataframe pandas a partir de um arquivo sam"""
        check = self._input_check()
        if not check:
            raise WrongOrMissingInput(f"The file {self.filename} either could not "
                                      f"been found in input folder"
                                      f" or it's in the wrong format")
        self._soft_sam()
        path = f'{self.absolute_path}/input/soft_cliped.sam'

        with open(path, 'r') as f:
            lines = f.readlines()
        data = []
        for line in lines:
            fields = line.strip().split('\t')
            row = {
                'pos': int(fields[3]),
                'fpos': self._get_fpos(fields[5], int(fields[3]), fields[9]),
                'fmap': self._get_fmap(fields[5], int(fields[3])),
                'cigar': fields[5],
                'seq': fields[9]
            }
            data.append(row)

        df = pd.DataFrame(data)
        df_depth = pd.read_csv(f"{self.absolute_path}/input/depth.csv",
                               sep="\t", names=['read', 'position', 'depth'])
        os.remove(f"{self.absolute_path}/input/soft_cliped.sam")
        os.remove(f"{self.absolute_path}/input/depth.csv")
        return df, data, self.filename.replace('.bam', ''), df_depth

    def _input_check(self):
        """Valida o input"""
        path = f'{self.absolute_path}/input/{self.filename}'
        exist = os.path.exists(path)
        if exist and self.filename.endswith('.bam'):
            return True
        else:
            return False

    def _soft_sam(self):
        """Cria o arquivo soft_cliped.sam que possui apenas reads com softclip"""
        file = self.filename
        subprocess.call(['bash', f'{self.absolute_path}/input/parse.sh', file])

    @staticmethod
    def _get_fmap(cigar, pos):
        """Gera a posição final da região mapeada do read"""
        fmap = 0
        pattern = r'(\d+)(M)'
        soft_qt = re.findall(pattern, cigar)
        for tupla in soft_qt:
            fmap += int(tupla[0])

        fmap = fmap + pos

        return fmap

    def _get_fpos(self, cigar, pos, seq):
        first_soft, _, _ = self.get_soft_seq(cigar, seq)
        if type(first_soft) == str:
            new_len = len(seq) - len(first_soft)
            fpos = new_len + pos
            return fpos
        else:
            fpos = len(seq) + pos
            return fpos

    @staticmethod
    def get_soft_seq(cigar, seq):
        """Gera 'strings' com a região inicial e final de softclip
        e também com a última região mapeada"""
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


class OutputHandler:
    """Gera o arquivo de log. Processará futuros outputs"""
    def __init__(self, absolute_path):
        self.absolute_path = absolute_path

    def create_log(self, text: str, mode: str):
        with open(f'{self.absolute_path}/output/log.txt', f'{mode}') as f:
            f.write(text)


class WrongOrMissingInput(Exception):
    pass
