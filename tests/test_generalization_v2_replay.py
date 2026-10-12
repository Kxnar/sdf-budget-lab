import sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import numpy as np
from generalization_v2_reproduce_amended3 import baseline_ray_evidence
class ReplayCheck(unittest.TestCase):
    def test_complete_events_and_signed_witnesses(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);ref=folder/'ref';ref.mkdir();cap=8192
            (folder/'fields.csv').write_text('cap,grid_n\n8192,3\n',encoding='utf-8')
            (ref/'rays.txt').write_text('0 0 0 0 1 0 0 -1 1\n0 0 0 0 1 0 0 -1 1\n',encoding='utf-8')
            (ref/'reference_roots.csv').write_text('ray,roots\n0,-0.25;0.25\n1,-0.25;0.25\n',encoding='utf-8')
            values=np.broadcast_to(np.array([.75,-.25,.75],dtype='<f4'),(3,3,3)).copy();values.tofile(folder/f'{cap}.bin')
            events=baseline_ray_evidence(folder,ref,cap,'grid')[0]
            np.testing.assert_allclose(events['roots'],[-.25,.25]);self.assertEqual(events['tangencies'],[]);self.assertFalse(events['ambiguous'])
            self.assertTrue(events['signed_witnesses'][1]['values'][3]<0)
            values[1,1,1]=-.5;values.tofile(folder/f'{cap}.bin')
            self.assertNotEqual(events,baseline_ray_evidence(folder,ref,cap,'grid')[0])
if __name__=='__main__':unittest.main()
