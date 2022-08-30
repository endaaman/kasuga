from collections import OrderedDict, namedtuple
import itertools

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import seaborn as sns

from endaaman import Commander


class C(Commander):
    def arg_hi(self, parser):
        # parser.add_argument('--target', '-t', default='br')
        pass

    def run_hi(self):
        print('hi')


c = C()
c.run()
