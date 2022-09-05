from collections import OrderedDict, namedtuple
import itertools
from itertools import combinations

import numpy as np
import pandas as pd
import pickle
import seaborn as sns
from matplotlib import pyplot as plt
from sklearn.linear_model import LinearRegression
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
import statsmodels.api as sma

from endaaman import Commander


class C(Commander):
    def run_cache(self):
        df = pd.read_csv('./data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct', sep='\t', skiprows=2, index_col=0)
        # print(df.shape)
        # In [9]: np.where(df.index.isin(['FBXO11',  'CIITA', 'HLA-DRA', ]))
        # Out[9]: (array([ 5820, 18006, 41214]),)

        # retrieve only 'FBXO11',  'CIITA', 'HLA-DRA'
        skip = np.arange(56203)
        # skip = np.delete(skip, [2, 5822, 18008, 41216])
        skip = np.delete(skip, [2, 5823, 18009, 41217])
        df = pd.read_csv('./data/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct', sep='\t', skiprows=skip, index_col=1)
        df.to_pickle('data/gtex/v8cache.pickle')

    def pre_common(self):
        with open('data/gtex/v8cache.pickle', 'rb') as f:
            df = pickle.load(f)
        df = df.drop(columns=['Name']).T

        self.ll = pd.read_csv('./data/gtex/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', sep='\t',  index_col=0)['SMTS']

        mm = pd.merge(df, self.ll, how='left', left_index=True, right_index=True)

        self.target_cols = ['FBXO11', 'CIITA', 'HLA-DRA']
        tmps = []
        for col in self.target_cols:
            tmp = mm[['SMTS', col]]
            tmp = tmp.rename(columns={col: 'value'})
            tmp['gene'] = col
            tmps.append(tmp)

        self.vv = pd.concat(tmps)
        self.vv['value'] = np.log2(self.vv['value'] )


    def run_violin(self):
        ''' 3遺伝子セットで組織ごとにプロット
        '''
        fig = plt.figure()
        split_count = 4
        for i, cols in enumerate(np.array_split(vv['SMTS'].unique(), split_count)):
            ax = fig.add_subplot(split_count, 1, i + 1)
            v = vv[vv['SMTS'].isin(cols)]
            sns.violinplot(x='SMTS', y='value', data=v, hue='gene', dodge=True,
                          jitter=True, color='black', palette='Set3', ax=ax)
            ax.set_ylabel('log2(TPM)')
            ax.tick_params(axis='x', labelrotation=45)
            ax.legend(loc='center left', bbox_to_anchor=(1, 0.7))

        plt.subplots_adjust(hspace=0.6)
        plt.show()

    def run_violin(self):
        ''' 遺伝子別に組織ごとにプロット
        '''
        fig = plt.figure()
        for i, col in enumerate(self.target_cols):
            ax = fig.add_subplot(len(self.target_cols), 1, i + 1)
            v = vv[vv['gene'] == col]
            ax.set_title(col)
            sns.violinplot(x='SMTS', y='value', data=v, dodge=True,
                          jitter=True, color='black', palette='Set3', ax=ax)
            ax.set_ylabel('log2(TPM)')
            ax.tick_params(axis='x', labelrotation=45)

        plt.subplots_adjust(hspace=0.6)
        plt.show()


    def run_coef(self):
        '''部位ごとの3遺伝子の相関、p valueを算出
        '''
        ii = self.vv.groupby('SMTS')

        data = []
        for (i, v) in ii:
            row = OrderedDict({
                'name': i,
            })
            for (a, b) in combinations(self.target_cols, 2):
                aa = v[v['gene'] == a]['value'].values
                bb = v[v['gene'] == b]['value'].values

                lr = LinearRegression()
                lr.fit(aa.reshape(-1, 1), bb)
                # row[f'{a} vs {b} coef'] = f'{lr.coef_[0]:.3f}'

                # X = aa.reshape(-1, 1)
                # y = bb
                # X_ = np.append(np.ones((len(X),1)), X, axis=1)
                # theta = np.append(lr.intercept_, lr.coef_)
                # y_preds = lr.predict(X)
                # RSS = np.sum((y-y_preds)**2)
                # RSE = np.sqrt(RSS/(len(X_)-len(X_[0])))
                # SE_sq = RSE**2 * np.linalg.inv(np.dot(X_.T,X_)).diagonal()
                # t = theta/np.sqrt(SE_sq)
                # p = [2*(1-stats.t.cdf(np.abs(t_val),(len(X_)-len(X_[0])))) for t_val in t]
                # row[f'{a} vs {b} p-value2'] = f'{p[1]:.6f}'

                aa2 = sma.add_constant(aa)
                est = sma.OLS(bb, aa2)
                r = est.fit()
                row[f'{a} vs {b} coef'] = f'{float(r.params[1]):.6f}'
                row[f'{a} vs {b} p-value'] = f'{float(r.pvalues[1]):.6f}'

            data.append(row)

        df = pd.DataFrame(data)
        # p = 'out/coef.xlsx'
        # df.to_excel(p)
        p = 'out/coef.tsv'
        df.to_csv(p, sep='\t')
        print(f'saved {p}')

if __name__ == '__main__':
    c = C()
    c.run()
