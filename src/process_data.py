import subprocess
import pandas as pd
import os


def create_data(filename):
    """Cria um dataframe pandas a partir de um arquivo sam"""
    check = input_check(filename)
    if check:
        soft_sam(filename)
        path = 'input/soft_cliped.sam'

        with open(path, 'r') as f:
            lines = f.readlines()
        data = []
        for line in lines:
            fields = line.strip().split('\t')
            row = {
                'qname': fields[0],
                'flag': int(fields[1]),
                'rname': fields[2],
                'pos': int(fields[3]),
                'fpos': int(fields[3]) + len(fields[9]),
                'cigar': fields[5],
                'seq': fields[9],
                'qual': fields[10],
            }
            data.append(row)

        df = pd.DataFrame(data)
        os.remove("input/soft_cliped.sam")
        return df
    else:
        raise WrongOrMissingInput(f"The file {filename} either could not "
                                  f"been found in input folder"
                                  f" or it's in the wrong format")


def input_check(filename):
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


class WrongOrMissingInput(Exception):
    pass
