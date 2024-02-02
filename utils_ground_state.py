
# find ground state and write to file
def find_ground_state(my_eigenvectors,my_min_index,append_mode,write_eigenstate_file=0):
    psi_g=my_eigenvectors[:,my_min_index]

    array_to_latex(psi_g,max_size=16)

    psi_g
    Xarr=np.zeros((mmax,2) )

    for i in range(mmax):
        Xarr[i][0]=np.real(psi_g[i])
        Xarr[i][1]=np.imag(psi_g[i])


    if write_eigenstate_file==1:
        if append_mode==0:
            filename_eigenstate="f2/eigenstate_data/f2.exact.L"+str(L)+"_state.msqr"+str(msqr)+".Gsqr"+str(Gsqr)
            df=pd.DataFrame(Xarr,dtype=np.float64) 
            #df.to_csv(filename_eigenstate,index=False,header=0,mode='w')
            append_mode=1
        
        else:
            filename_eigenstate="f2/eigenstate_data/f2.exact.L"+str(L)+"_state.msqr"+str(msqr)+".Gsqr"+str(Gsqr)
            df=pd.DataFrame(Xarr,dtype=np.float64) 
            #df.to_csv(filename_eigenstate,index=False,header=0,mode='a')

        Yarr=np.array(pd.read_csv(filename_eigenstate))
        Yarr
        return psi_g,Yarr
    else:
        return psi_g,[]