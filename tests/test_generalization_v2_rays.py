import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from generalization_v2_rays import *
class SegmentedRayTests(unittest.TestCase):
    def test_narrow_pair(self):
        roots,_,_=polynomial_roots([.25-1e-12,-1,1,0]);self.assertEqual(len(roots),2);self.assertLess(roots[1]-roots[0],3e-6)
    def test_tangent(self):
        roots,tangencies,_=polynomial_roots([.25,-1,1,0]);self.assertEqual(roots,[]);self.assertEqual(tangencies,[.5])
    def test_three_crossings(self):
        roots,_,_=polynomial_roots([-.09375,.6875,-1.5,1]);np.testing.assert_allclose(roots,[.25,.5,.75],atol=3e-15)
    def test_near_tangent_signs(self):
        self.assertEqual(len(polynomial_roots([.25+1e-16,-1,1,0])[0]),0);self.assertEqual(len(polynomial_roots([.25-1e-16,-1,1,0])[0]),2)
    def test_nearly_repeated_cubic(self):
        c=np.polynomial.polynomial.polymul(np.polynomial.polynomial.polymul([-.5,1],[-.500001,1]),[-.8,1]);self.assertEqual(len(polynomial_roots(c)[0]),3)
    def test_zero(self):self.assertTrue(polynomial_roots([0,0,0,0])[2])
    def test_endpoints_not_interior(self):self.assertEqual(polynomial_roots([0,-1,1,0])[0],[])
    def test_reversed_signs_fail(self):
        ray=np.array([0,0,0,0,1,0,0,-1,1]);self.assertFalse(signed_intervals(lambda p:np.where(abs(p[:,0])<.5,1.,-1.),ray,np.array([-.5,.5]),np.array([-.5,.5])))
    def test_decoded_precision(self):
        values=np.array([[.1,.2,.3,.4,.5,.6,.7,.8]]);q=decoded_corners(values,'q8');self.assertTrue(np.all(np.isfinite(q)));self.assertTrue(np.array_equal(decoded_corners(values,'f32'),values.astype(np.float32).astype(float)))
    def test_higher_order_interface_sides(self):
        ray=np.array([0,0,0,0,1,1,1,0,1]);corners=np.array([[bool(c&1),bool(c&2),bool(c&4)] for c in range(8)],dtype=float)
        for degree in (2,3):
            values=np.prod(corners[:,:degree]-.5,axis=1);segments=np.array([[0,.5,0,0,0,1,1,1,*values],[.5,1,0,0,0,1,1,1,*values]])
            result=transitions(segments,ray)
            if degree==2:self.assertEqual(len(result['roots']),0);np.testing.assert_array_equal(result['tangencies'],[.5])
            else:np.testing.assert_array_equal(result['roots'],[.5]);self.assertEqual(len(result['jumps']),0)
            self.assertFalse(result['ambiguous'])
if __name__=='__main__':unittest.main()
