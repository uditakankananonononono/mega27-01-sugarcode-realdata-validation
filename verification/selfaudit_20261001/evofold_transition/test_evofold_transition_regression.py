import numpy as np
from sugarcode.modules.evofold_4d import transition_trace
C=np.array([[0,0,0],[2,0,0],[0,3,0],[0,0,4],[2,3,1],[1,2,4]],float)
def reference():
 n=len(C); H=np.zeros((3*n,3*n))
 for i in range(n):
  for j in range(i+1,n):
   d=C[j]-C[i];b=-np.outer(d,d)/(d@d)
   H[3*i:3*i+3,3*j:3*j+3]=b;H[3*j:3*j+3,3*i:3*i+3]=b
   H[3*i:3*i+3,3*i:3*i+3]-=b;H[3*j:3*j+3,3*j:3*j+3]-=b
 w,V=np.linalg.eigh(H);return H,w,V[:,np.where(w>1e-6)[0][0]]
def test_reference_has_six_rigid_body_modes():
 H,w,v=reference(); assert sum(abs(w)<1e-6)==6
 assert np.max(abs(H@np.tile([1.,0,0],len(C))))<1e-12
 assert np.linalg.norm(v.reshape(-1,3).mean(0))<1e-12
def test_transition_nonrigid_mode_preserves_centroid():
 r=transition_trace(C,steps=3,amplitude=3);d=np.array(r['frames'][1])-C
 assert np.linalg.norm(d.mean(0))<.002

def test_transition_uses_first_nonrigid_anm_mode_up_to_sign():
 H,w,v=reference();r=transition_trace(C,steps=3,amplitude=3)
 d=(np.array(r['frames'][1])-C).ravel();d/=np.linalg.norm(d)
 assert abs(d@v)>.999
