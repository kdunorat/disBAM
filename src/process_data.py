import subprocess
import pandas as pd
import os
import re


def create_data(filename):
    """Cria um dataframe pandas a partir de um arquivo sam"""
    check = input_check(filename)
    if not check:
        raise WrongOrMissingInput(f"The file {filename} either could not "
                                  f"been found in input folder"
                                  f" or it's in the wrong format")
    soft_sam(filename)
    path = 'input/soft_cliped.sam'

    with open(path, 'r') as f:
        lines = f.readlines()
    data = []
    for line in lines:
        fields = line.strip().split('\t')
        row = {
            'pos': int(fields[3]),
            'fpos': get_fpos(fields[5], int(fields[3]), fields[9]),
            'fmap': get_fmap(fields[5], int(fields[3])),
            'cigar': fields[5],
            'seq': fields[9]
        }
        data.append(row)

    df = pd.DataFrame(data)
    os.remove("input/soft_cliped.sam")
    return df, data


def input_check(filename):
    """Valida o input"""
    path = f'{os.getcwd()}/input/{filename}'
    exist = os.path.exists(path)
    if exist and filename.endswith('.bam'):
        return True
    else:
        return False


def soft_sam(filename):
    """Cria o arquivo soft_cliped.sam que possui apenas reads com softclip"""
    file = filename
    subprocess.call(['bash', 'input/parse.sh', file])


def get_fmap(cigar, pos):
    """Gera a posição final da região mapeada do read"""
    fmap = 0
    pattern = r'(\d+)(M)'
    soft_qt = re.findall(pattern, cigar)
    for tupla in soft_qt:
        fmap += int(tupla[0])

    fmap = fmap + pos

    return fmap


def get_fpos(cigar, pos, seq):
    first_soft, _, _ = get_soft_seq(cigar, seq)
    if first_soft:
        new_len = len(seq) - len(first_soft)
        fpos = new_len + pos
        return fpos
    else:
        fpos = len(seq) + pos
        return fpos


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


class WrongOrMissingInput(Exception):
    pass
