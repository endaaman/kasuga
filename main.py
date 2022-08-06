from collections import OrderedDict
import itertools

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
# import seaborn as sns

from endaaman import Commander

db_names = OrderedDict({
    'brca': 'data/TCGA.BRCA.sampleMap_HiSeqV2.txt',
    'brca_exon': 'data/TCGA.BRCA.sampleMap_HiSeqV2_exon.txt',
    'brca_pancan': 'data/TCGA.BRCA.sampleMap_HiSeqV2_PANCAN.txt',
    'paad': 'data/TCGA.PAAD.sampleMap_HiSeqV2.txt',
    'paad_pancan': 'data/TCGA.PAAD.sampleMap_HiSeqV2_PANCAN.txt',
})

required_colums = ['FBXO11', 'HLA-DRA', 'CIITA']
target_colums = ['FBXO11', 'HLA-DRA', 'CIITA']

def load_hiseq_data(name):
    df = pd.read_csv(db_names[name], sep='\t', index_col=0,).T
    print(f'loaded {name}')
    return df


class C(Commander):

    def run_hi(self):
        df = load_hiseq_data('brca')

        pairs = list(itertools.combinations(target_colums, 2))
        fig = plt.figure(figsize=(8, 16))

        for i, (a, b) in enumerate(pairs):
            ax = fig.add_subplot(len(pairs), 1, i + 1)
            df_a, df_b  = df[a], df[b]
            r = np.corrcoef(df_a, df_b)[0][1]
            ax.set_title(f'{a} vs {b} coef={r:.2f}')
            ax.scatter(df_a, df_b)

        plt.show()

    def run_3d(self):
        df = load_hiseq_data('brca')

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        ax.scatter3D(df['FBXO11'], df['CIITA'], df['HLA-DRA'])
        ax.set_xlabel('FBXO11')
        ax.set_ylabel('CIITA')
        ax.set_zlabel('HLA-DRA')

        plt.show()



c = C()
c.run()
