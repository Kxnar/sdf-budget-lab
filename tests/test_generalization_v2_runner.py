import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from generalization_v2_run import expected_jobs,FIT_SEEDS,COMPACT,GUARD
from generalization_common import METHODS,CAPS
from generalization_v2_report import crossing

class FinalMatrixTests(unittest.TestCase):
    def test_cartesian_deterministic_reuse(self):
        f=dict(geometries=[dict(id=str(i)) for i in range(4)],methods=METHODS+list(COMPACT),fit_seeds=FIT_SEEDS)
        jobs=expected_jobs(f)
        self.assertEqual(len(jobs),156);self.assertEqual(len(jobs)*len(CAPS),2028)
        self.assertEqual(sum(m in COMPACT for g,s,m in jobs),12)
        self.assertEqual(sum(m=='grid' for g,s,m in jobs),4)
        self.assertEqual(sum(m not in COMPACT and m!='grid' for g,s,m in jobs),140)
        self.assertEqual(GUARD,1610612736)
    def test_nonmonotonic_actual_crossing_and_censor(self):
        rows=[dict(geometry='g',method='m',seed=340,cap=c,bytes=b,**{'pass':p}) for c,b,p in [(8192,8000,False),(12288,11000,True),(16384,10000,True),(24576,24000,False)]]
        self.assertEqual(crossing(rows,'g','m',340)['bytes'],10000)
        self.assertIsNone(crossing(rows,'missing','m',340))
if __name__=='__main__':unittest.main()
