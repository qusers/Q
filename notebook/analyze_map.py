import sys
import os
from pathlib import Path

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/')))

import QmapFEP

a = QmapFEP.Analyze('/home/willem/QmapFEP-test/test.json',
"/home/willem/software/Q/data/benchmark/qligfep/raw/1.JACS/OPLS2015/2.CDK2/no-vs/2fs")