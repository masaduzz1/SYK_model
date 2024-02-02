# Initial state preperation
def edge_flipped_instate(Q):
    myzero=np.array([1,0])
    myone=np.array([0,1])
    
    psi=myone
    
    for i in range(1,Q-1):
        psi=np.kron(psi,myzero)
        
    
    psi=np.kron(psi,myone)
    
    psi= psi +0.0j
    
    return psi
    
def left_edge_flipped_instate(Q):
    myzero=np.array([1,0])
    myone=np.array([0,1])
    
    psi=myone
    
    for i in range(1,Q):
        psi=np.kron(psi,myzero)
        
    
    
    psi= psi +0.0j
    
    return psi
    
left_edge_flipped_instate(3) # 01 10 01  --> 01 0100 -> 0000_0100
def center_flipped_instate(Q):
    myzero=np.array([1,0])
    myone=np.array([0,1])
    
    psi=myzero
    
    for i in range(1,Q):
        if i==Q//2:
            psi=np.kron(psi,myone)
        else:
            psi=np.kron(psi,myzero)
        
    
    
    psi= psi +0.0j
    
    return psi
def all_zero_qubit(Q):
    psi=np.zeros(2**Q)+0.0j
    psi[0]=1
    return psi



def random_unitary_create(Q):
    u_op= unitary_group.rvs(2**Q)
    instate=all_zero_qubit(Q)
    #display(array_to_latex(u_op))
    instate=u_op @ instate
    #display(array_to_latex(instate))
    return instate

def my_ceu(ce,Q):
    mmax=2**Q
    u_op= ce.gen_cue(mmax)
    instate=all_zero_qubit(Q)
    #display(array_to_latex(u_op))
    instate=u_op @ instate
    #display(array_to_latex(instate))
    return instate

#ce=Circular()
#for _ in range(5):
#    array_to_latex(my_ceu(ce,2))