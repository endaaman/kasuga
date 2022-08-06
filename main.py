from collections import OrderedDict, namedtuple
import itertools

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
# import seaborn as sns

from endaaman import Commander

Collection = namedtuple('Collection', ['abbr', 'name', 'path'])

collections = OrderedDict({
    c.abbr:c for c in [
        Collection('br', 'Breast cancer', 'data/TCGA.BRCA.sampleMap_HiSeqV2.txt'),
        Collection('pan', 'Pancreatic cancer', 'data/TCGA.PAAD.sampleMap_HiSeqV2.txt'),
        Collection('colrec',  'Collon rectum cancer', 'data/TCGA.COADREAD.sampleMap_HiSeqV2.txt'),
    ]
})

required_colums = ['FBXO11', 'HLA-DRA', 'CIITA']
target_colums = ['FBXO11', 'HLA-DRA', 'CIITA']


def load_hiseq_data(c):
    df = pd.read_csv(c.path, sep='\t', index_col=0,).T
    print(f'loaded {c.name} data')
    return df


class C(Commander):
    def arg_hi(self, parser):
        parser.add_argument('--target', '-t', default='br')

    def run_hi(self):
        c = collections[self.args.target]
        df = load_hiseq_data(c)

        pairs = list(itertools.combinations(target_colums, 2))
        fig = plt.figure(figsize=(8, 14))

        for i, (a, b) in enumerate(pairs):
            ax = fig.add_subplot(len(pairs), 1, i + 1)
            df_a, df_b  = df[a], df[b]
            r = np.corrcoef(df_a, df_b)[0][1]
            ax.set_title(f'{a} vs {b} coef={r:.2f}')
            ax.scatter(df_a, df_b)

        plt.suptitle(f'{c.name}')
        plt.show()
        plt.savefig('out/{c.abbr}.png')

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
