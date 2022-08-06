from collections import OrderedDict

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

def load_hiseq_data(name):
    df = pd.read_csv(db_names[name], sep='\t', index_col=0,).T
    print(f'loaded {name}')
    return df


class C(Commander):

    def run_hi(self):
        df = load_hiseq_data('brca')

        n = 3
        fig = plt.figure(figsize=(8, 16))
        ax = fig.add_subplot(n, 1, 1)
        ax.set_title('FBXO11 vs HLA-DRA')
        ax.scatter(df['FBXO11'], df['HLA-DRA'])

        ax = fig.add_subplot(n, 1, 2)
        ax.set_title('FBXO11 vs CIITA')
        ax.scatter(df['FBXO11'], df['CIITA'])

        ax = fig.add_subplot(n, 1, 3)
        ax.set_title('HLA-DRA vs CIITA')
        ax.scatter(df['HLA-DRA'], df['CIITA'])
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
