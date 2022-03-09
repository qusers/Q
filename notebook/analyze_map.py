import sys
import os

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/')))

import QmapFEP

a = QmapFEP.Analyze('/home/willem/QmapFEP-test/test.json',
"/home/willem/software/Q/data/benchmark/qligfep/raw/1.JACS/OPLS2015/2.CDK2/no-vs/2fs")

b = QmapFEP.LoadData('/home/willem/QmapFEP-test/test.json','/home/willem/software/Q/data/benchmark/qligfep/processed/CDK2_results.csv')
c = QmapFEP.GenPlot('/home/willem/QmapFEP-test/test.json','/home/willem/QmapFEP-test')