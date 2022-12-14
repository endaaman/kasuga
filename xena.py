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
        Collection('brca', 'Breast cancer', 'data/xena/TCGA.BRCA.sampleMap_HiSeqV2.txt'),
        Collection('paad', 'Pancreatic cancer', 'data/xena/TCGA.PAAD.sampleMap_HiSeqV2.txt'),
        Collection('colrec',  'Colon rectum cancer', 'data/xena/TCGA.COADREAD.sampleMap_HiSeqV2.txt'),
        Collection('gbm',  'Glioblastoma', 'data/xena/TCGA.COADREAD.sampleMap_HiSeqV2.txt'),
    ]
})

required_colums = ['FBXO11', 'HLA-DRA', 'CIITA']
target_colums = ['FBXO11', 'HLA-DRA', 'CIITA']


def load_xena_data(name):
    c = collections[name]
    df = pd.read_csv(c.path, sep='\t', index_col=0).T
    print(f'loaded {c.name} data')
    return df


class C(Commander):
    def arg_hi(self, parser):
        parser.add_argument('--target', '-t', default='br')

    def run_hi(self):
        df = load_xena_data(self.args.target)

        pairs = list(itertools.combinations(target_colums, 2))
        fig = plt.figure(figsize=(8, 14))

        for i, (a, b) in enumerate(pairs):
            ax = fig.add_subplot(len(pairs), 1, i + 1)
            df_a, df_b  = df[a], df[b]
            r = np.corrcoef(df_a, df_b)[0][1]
            ax.set_title(f'{a} vs {b} coef={r:.2f}')
            ax.scatter(df_a, df_b)

        plt.suptitle(f'{c.name} n={df.shape[0]}')
        plt.show()
        plt.savefig(f'out/{c.abbr}.png')

    def run_3d(self):
        df = load_xena_data('brca')

        fig = plt.figure(figsize=(8, 8))
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        ax.scatter3D(df['FBXO11'], df['CIITA'], df['HLA-DRA'])
        ax.set_xlabel('FBXO11')
        ax.set_ylabel('CIITA')
        ax.set_zlabel('HLA-DRA')

        plt.show()


if __name__ == '__main__':
    c = C()
    c.run()
