from scipy import linalg as sLA
import numpy as np
from numpy import linalg as LA



def density_matrix(myH,mybeta):
    rho=sLA.expm(-mybeta*myH)
    return rho/( np.real( np.trace(rho)) )

def time_evolve(tt,H):
    return sLA.expm(-1.0j*tt*H)
    
def heisenberg_op(Ut,W): 
    return (np.conj(Ut).T) @ W @ Ut 


#compute Tr(rho^1/4 W rho^1/4 V  rho^1/4 W rho^1/4 V)
def otoc(rho_beta4,tt,syk_op,W,V):
    #display(rho_beta4)
    
    Ut=time_evolve(tt,syk_op)
    #display(array_to_latex(Ut))
    
    W_t=heisenberg_op(Ut,W) 
    #display(array_to_latex(W_t) )
    
    
    rho_W_t=rho_beta4 @ W_t
    rho_Vop=rho_beta4 @  V

    RWt_RV_RWt_RV= LA.matrix_power( (rho_W_t @ rho_Vop ),2)  
    
    #RWt_RVdag_RWt_RV=(rho_W_t @ rho_beta4 @ np.conj(V.T) ) @ (rho_W_t @ rho_Vop )

    A=np.real(np.trace(RWt_RV_RWt_RV))

    #denominator 
    #W_t_sqr_V_sqr= (rho_W_t @ rho_W_t)  @ ( rho_Vop @ rho_Vop)
    #print("normalization")
    #array_to_latex(W_t_sqr_Vdag_V,max_size=16)
    #B=np.real(np.trace( W_t_sqr_V_sqr))
    #display(array_to_latex(RWt_RV_RWt_RV))
    return A

