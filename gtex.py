import os
import pickle
from collections import OrderedDict, namedtuple
from itertools import combinations

import numpy as np
import pandas as pd
from tqdm import tqdm
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
import statsmodels.api as sma

from endaaman import Commander


class C(Commander):
    def make_v8_kasuga_cache(self):
        # df = pd.read_csv('./data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct', sep='\t', skiprows=2, index_col=0)
        # 56200 genes
        # print(df.shape)
        # In [9]: np.where(df.index.isin(['FBXO11',  'CIITA', 'HLA-DRA', ]))
        # Out[9]: (array([ 5820, 18006, 41214]),)

        # retrieve only 'FBXO11',  'CIITA', 'HLA-DRA'
        skip = np.arange(56203)
        # skip = np.delete(skip, [2, 5822, 18008, 41216])
        skip = np.delete(skip, [2, 5823, 18009, 41217])
        print('loading original data..')
        df = pd.read_csv('./data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct', sep='\t', skiprows=skip, index_col=1)
        df.to_pickle('data/gtex/v8_kasuga_cache.pickle')
        return df

    target_cols = ['FBXO11', 'CIITA', 'HLA-DRA']

    def load_v8_kasuga(self):
        p = 'data/gtex/v8_kasuga_cache.pickle'
        if os.path.exists(p):
            df = pickle.load(open(p, 'rb'))
        else:
            df = self.make_v8_kasuga_cache()

        values = df.drop(columns=['Name']).T
        labels = pd.read_csv('./data/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', sep='\t',  index_col=0)[['SMTS']]

        merged = pd.merge(values, labels['SMTS'], how='left', left_index=True, right_index=True)
        tmps = []
        for col in self.target_cols:
            tmp = merged[['SMTS', col]]
            tmp = tmp.rename(columns={col: 'value'})
            tmp['gene'] = col
            tmps.append(tmp)

        vv = pd.concat(tmps)
        vv['value'] = np.log2(vv['value'])

        return values, labels, vv


    def arg_violin(self, parser):
        parser.add_argument('-b', '--by', choices=['tissue', 'gene'], default='tissue')

    def run_violin(self):
        ''' 3遺伝子セットで組織ごとにプロット
        '''
        _, _, vv = self.load_v8_kasuga()
        vv.to_csv('out/violin.tsv', sep='\t')
        fig = plt.figure(figsize=(20, 10))

        if self.args.by == 'tissue':
            split_count = 4
            for i, cols in enumerate(np.array_split(vv['SMTS'].unique(), split_count)):
                ax = fig.add_subplot(split_count, 1, i + 1)
                v = vv[vv['SMTS'].isin(cols)]
                sns.violinplot(x='SMTS', y='value', data=v, hue='gene', dodge=True,
                              jitter=True, color='black', palette='Set3', ax=ax)
                ax.set_ylabel('log2(TPM)')
                ax.tick_params(axis='x', labelrotation=20)
                ax.set_xlabel('')
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.7))

        elif self.args.by == 'gene':
            for i, col in enumerate(self.target_cols):
                ax = fig.add_subplot(len(self.target_cols), 1, i + 1)
                v = vv[vv['gene'] == col]
                ax.set_title(col)
                sns.violinplot(x='SMTS', y='value', data=v, dodge=True,
                              jitter=True, color='black', palette='Set3', ax=ax)
                ax.set_ylabel('log2(TPM)')
                ax.tick_params(axis='x', labelrotation=45)
                ax.set_xlabel('')

        plt.subplots_adjust(hspace=0.6)
        plt.savefig(f'out/violin_by_{self.args.by}.png')
        plt.show()

    def run_coef(self):
        '''部位ごとの3遺伝子の相関、p valueを算出
        '''
        _, _, vv = self.load_v8_kasuga()
        ii = vv.groupby('SMTS')

        data = []
        for (i, v) in ii:
            row = OrderedDict({
                'name': i,
            })
            for (a, b) in combinations(self.target_cols, 2):
                aa = v[v['gene'] == a]['value']
                bb = v[v['gene'] == b]['value']

                lr = LinearRegression()
                lr.fit(aa.values.reshape(-1, 1), bb.values)
                row[f'{a} vs {b} coef'] = f'{lr.coef_[0]:.3f}'

                corr = aa.corr(bb)
                row[f'{a} vs {b} corr'] = f'{corr:.3f}'


            data.append(row)

        df = pd.DataFrame(data)
        # p = 'out/coef.xlsx'
        # df.to_excel(p)
        p = 'out/coef.tsv'
        df.to_csv(p, sep='\t')
        print(f'saved {p}')

    def run_calc_coef(self):
        df = pd.read_csv('data/gtex/gene_reads_2017-06-05_v8_liver.gct', sep='\t', skiprows=2, index_col=2)
        df.drop(columns=['id', 'Name'], inplace=True)
        self.df = df
        return

        t = tqdm(list(combinations(df.columns, 2)))
        for (a, b) in t:
            pass
            # t.set_description(f'{a} {b}')
            # t.refresh()
            # print(a, b)
            # df[a]
            # df[b]

if __name__ == '__main__':
    c = C()
    c.run()
