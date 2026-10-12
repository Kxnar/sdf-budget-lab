"""Meaningful gate, ray, byte-layout and resume/censoring regression tests."""
import os,sys,tempfile,unittest
import subprocess
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from generalization_common import *
from generalization_verify import strict_gate,ray_sign_pattern
from generalization_run import expected_jobs
from distortion_analysis import layout,decode
class HarnessTests(unittest.TestCase):
    def test_and_gate(self):
        self.assertTrue(strict_gate(True,True,True,True,True))
        for i in range(5):
            flags=[True]*5;flags[i]=False;self.assertFalse(strict_gate(*flags))
    def test_width_missing_extra_discontinuous(self):
        rays=np.array([[0,0,0,0,1,0,0,-1,1]])
        roots,continuous=independent_roots(lambda p:np.abs(p[:,0])-.1,rays)[0];np.testing.assert_allclose(roots,[-.1,.1],atol=1e-10);self.assertTrue(continuous)
        missing,_=independent_roots(lambda p:np.ones(len(p)),rays)[0];self.assertEqual(len(missing),0)
        extra,_=independent_roots(lambda p:np.minimum(np.abs(p[:,0]-.3),np.abs(p[:,0]+.3))-.1,rays)[0];self.assertEqual(len(extra),4)
        roots,continuous=independent_roots(lambda p:np.where(np.abs(p[:,0])<.1,-1.,1.),rays)[0];self.assertFalse(continuous);np.testing.assert_allclose(roots,[-.1,.1],atol=1e-10)
    def test_reversed_interval_signs_with_identical_roots_and_widths(self):
        ray=np.array([0,0,0,0,1,0,0,-1,1]);roots=np.array([-.1,.1]);correct=lambda p:np.abs(p[:,0])-.1;reversed_sign=lambda p:-correct(p)
        recovered,_=independent_roots(reversed_sign,ray[None,:])[0];np.testing.assert_allclose(recovered,roots,atol=1e-10);np.testing.assert_allclose(np.diff(recovered),np.diff(roots),atol=1e-10);self.assertTrue(ray_sign_pattern(correct,ray,roots,roots));self.assertFalse(ray_sign_pattern(reversed_sign,ray,roots,recovered))
    def test_nonmonotonic_actual_crossing_and_censor(self):
        rr=[{'bytes':70,'pass':False},{'bytes':90,'pass':True},{'bytes':80,'pass':True},{'bytes':100,'pass':False}];self.assertEqual(crossing(rr),80);self.assertIsNone(crossing([{'bytes':262144,'pass':False}]))
    def test_hash_reject_changed_resume(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'field.bin';path.write_bytes(b'original');m=manifest(Path(d));check_hashes(Path(d),m);path.write_bytes(b'changed')
            with self.assertRaises(AssertionError):check_hashes(Path(d),m)
    def test_exact_layout_and_overlap_rejection(self):
        # A root q8 leaf: 4-byte descriptor +16-byte payload +64-byte object.
        raw=struct.pack('<Iff8B',7,-.1,.001,*range(8));layout(raw,'mse',5);self.assertEqual(len(raw)+64,84);p=decode(raw,np.array([[0,0,0]]));self.assertTrue(np.isfinite(p).all())
        with self.assertRaises(AssertionError):layout(struct.pack('<Iff8B',3,-.1,.001,*range(8)),'mse',5)
    def test_cartesian_completeness(self):
        jobs=expected_jobs({'geometries':[{'id':'a'},{'id':'b'}]});self.assertEqual(len(jobs),72);keys={(g['id'],s,m) for g,s,m in jobs};self.assertEqual(len(keys),72);self.assertEqual(sum(m=='grid' for g,s,m in jobs),2)
    def test_depth_seven_decoder_terminal_and_layout(self):
        children={};nodes=1;current=0
        for depth in range(7):children[current]=nodes;current=nodes;nodes+=8
        raw=bytearray(nodes*4)
        for node in range(nodes):
            if node in children:struct.pack_into('<I',raw,node*4,children[node]*4)
            else:
                offset=len(raw);raw.extend(struct.pack('<ff8B',-.05 if node==current else 1.,0.,*[0]*8));struct.pack_into('<I',raw,node*4,offset|3)
        layout(bytes(raw),'candidate',7);values,_,depth=decode_packed(bytes(raw),np.array([[-.999,-.999,-.999],[.5,.5,.5]]),True);np.testing.assert_array_equal(depth,[7,1]);np.testing.assert_allclose(values,[-.05,1],atol=1e-8)
    def test_native_cache_identity_tail_and_conflict_rejection(self):
        exe=Path(os.environ.get('SDF_GENERALIZATION_BENCH_EXE',str(ROOT/'build-generalization-core/Release/sdf_generalization_bench.exe'))).resolve()
        if not exe.exists():
            if 'SDF_GENERALIZATION_BENCH_EXE' in os.environ:self.fail('Configured native test executable not available: '+str(exe))
            self.skipTest('Native test build not available')
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);spec=d/'geometry.txt';spec.write_text('test torus .035 0 0 .0175 1 0 0 0 1 0 0 0 1 0 0 0 ""\n');cache=d/'cache.bin';identity='a'*64
            for content,expected in [(identity.encode()+b'\ntruncated','Truncated shared cache'),(b'b'*64+b'\n','identity mismatch'),(identity.encode()+b'\n'+struct.pack('<4d',2,0,0,.1),'Invalid shared cache sample'),(identity.encode()+b'\n'+struct.pack('<4d',0,0,0,.1)+struct.pack('<4d',0,0,0,.2),'Contradictory shared cache labels')]:
                cache.write_bytes(content);result=subprocess.run([str(exe),'qualify',str(spec),str(d/'out'),'1','9','--cache',str(cache),'--cache-key',identity],capture_output=True,text=True);self.assertNotEqual(result.returncode,0);self.assertIn(expected,result.stderr)
if __name__=='__main__':unittest.main()
