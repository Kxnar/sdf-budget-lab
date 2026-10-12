import sys,tempfile,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from generalization_v2_reproduce_amended2 import same_acceptance
from generalization_common import check_hashes,sha
from generalization_v2_report_amended4 import failure_counts

class ReproductionChecks(unittest.TestCase):
    def test_only_timings_may_differ(self):
        same_acceptance([dict(roots='1;2',closed='1',seconds='3')],[dict(roots='1;2',closed='1',seconds='99')])
        for key,value in [('roots','1;3'),('closed','0'),('resolved','0')]:
            base=dict(roots='1;2',closed='1',resolved='1',seconds='3');changed=dict(base);changed[key]=value
            with self.assertRaises(AssertionError):same_acceptance([base],[changed])
    def test_stale_manifest_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            folder=Path(temp);p=folder/'flags.bin';p.write_bytes(b'\x01\x03');manifest={'flags.bin':sha(p)};check_hashes(folder,manifest);p.write_bytes(b'\x01\x01')
            with self.assertRaises(AssertionError):check_hashes(folder,manifest)
    def test_reference_unknown_does_not_hide_numerical(self):
        assessment=dict(status='reference_unresolved',scalar_pass=False,feature_pass=True,width_pass=False,metrics=[dict(unresolved=1,nonfinite=0)],topology=[])
        counts=failure_counts([{}]*13,assessment)
        self.assertEqual(counts['numerical_unresolved'],13)
        self.assertEqual(counts['scalar_observed_failed'],0)
    def test_decoded_invalid_gates_are_unassessed(self):
        counts=failure_counts([{}]*13,dict(status='decoded_invalid'))
        self.assertEqual(counts['scalar_feature_width_unassessed'],13)
        self.assertEqual(counts['width_observed_failed'],0)
        counts=failure_counts([{}],dict(scalar_pass=True,feature_pass=True,width_pass=True,metrics=[],topology=[dict(status='observed_failure')]))
        self.assertEqual(counts['observed_topology_failed'],1)
if __name__=='__main__':unittest.main()
