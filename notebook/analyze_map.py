import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/')))

import QmapFEP
b = QmapFEP.LoadData(o='../../tmp/QmapFEP-test/test.json',ifile='../data/benchmark/qligfep/processed/CDK2_results.csv')

a = QmapFEP.Analyze(o='../../tmp/QmapFEP-test/test.json',
datadir="../data/benchmark/qligfep/raw/1.JACS/OPLS2015/2.CDK2/no-vs/2fs", wd='../../tmp/QmapFEP-test')

c = QmapFEP.GenPlot(o='../../tmp/QmapFEP-test/test.json',wd='../../tmp/QmapFEP-test')
