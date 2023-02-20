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


class WrongOrMissingInput(Exception):
    pass
