"""Independent codec/layout and exact ray regression checks."""
import os,sys,tempfile,subprocess,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from compact_mesh import *

class Compact(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp=tempfile.TemporaryDirectory(prefix='compact-tests-')
        exe=Path(os.environ.get('SDF_COMPACT_TESTS_EXE',str(Path(__file__).resolve().parents[1]/'build-generalization-v2-core/Release/sdf_compact_tests.exe'))).resolve()
        subprocess.run([str(exe),cls.temp.name],check=True)
        cls.fields={i:decode((Path(cls.temp.name)/f'{i}.bin').read_bytes()) for i in range(1,6)}
    @classmethod
    def tearDownClass(cls):cls.temp.cleanup()
    def test_all_layouts_connectivity(self):
        vertices=np.array([[.5 if c&1 else -.5,.5 if c&2 else -.5,.5 if c&4 else -.5] for c in range(8)])
        faces=np.array([[0,2,3],[0,3,1],[4,5,7],[4,7,6],[0,1,5],[0,5,4],[2,6,7],[2,7,3],[0,4,6],[0,6,2],[1,3,7],[1,7,5]])
        for d in self.fields.values():self.assertEqual(source_check(d,vertices,faces),0);self.assertEqual(mesh_topology(d['vertices'],d['triangles']),(True,[2]))
    def test_exact_shared_event(self):
        ray=[0,0,0,0,1,0,0,-.9,.9]
        for d in self.fields.values():
            result=ray_events(d['vertices'],d['triangles'],ray);self.assertFalse(result['unresolved']);self.assertEqual(result['roots'],[-.5,.5])
    def test_tangent_and_disjoint_plane(self):
        d=self.fields[1]
        tangent=ray_events(d['vertices'],d['triangles'],[0,0,0,1,1,1,-1,.1,.9]);self.assertEqual(tangent['tangencies'],[.5]);self.assertEqual(tangent['roots'],[])
        disjoint=ray_events(d['vertices'],d['triangles'],[0,0,.5,0,1,0,0,.6,.9]);self.assertFalse(disjoint['unresolved'])
    def test_adjacent_crossings(self):
        d=self.fields[1];v=d['vertices'].copy();v[:,0]=np.where(v[:,0]<0,.5,np.nextafter(.5,1.))
        result=ray_events(v,d['triangles'],[0,0,0,0,1,0,0,.4,.9]);self.assertFalse(result['unresolved']);self.assertEqual(result['roots'],[.5,np.nextafter(.5,1.)])
    def test_xor_compressed(self):
        d=decode((Path(self.temp.name)/'xor-compressed.bin').read_bytes());self.assertTrue(any(b['mode']==1 for b in d['blocks']))
        v=np.array([[.4+(i%4)*2**-32,.4+((i//4)%4)*2**-32,.4+(i//16)*2**-32] for i in range(64)])
        self.assertEqual(source_check(d,v,np.array([[0,1,4]])),0)
    def test_truncated_all_codecs(self):
        for i in range(1,6):
            raw=(Path(self.temp.name)/f'{i}.bin').read_bytes()
            for length in (0,79,len(raw)-1):
                with self.assertRaises((AssertionError,struct.error)):decode(raw[:length])
if __name__=='__main__':unittest.main()
