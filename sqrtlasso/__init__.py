def sqrtlasso(X,Y,lam,options):
    MaxIter = options['MaxIter']
    PrintOut= options['PrintOut']
    OptTolNorm = options['OptTolNorm']
    
    n,p = X.shape

    centered_matrix = X - X.mean(axis=0)[np.newaxis,:]
    cov = np.dot(centered_matrix.T, centered_matrix)
    eigvals, eigvecs = np.linalg.eig(cov)
    gamma = np.dot(eigvecs, eigvals.T)
    gamma = abs(gamma)

    # start from ridge estimator
    RidgeMatrix = np.identity(p)

    for i in range(p):
        RidgeMatrix[j,j]= (lam)*gamma[j]

    a = (np.matmul(np.transpose(X),X) + RidgeMatrix)
    b = np.matmul(np.transpose(X),Y)
    beta = np.linalg.solve(a,b)


    # Start the iteration
    Iter  = 0

    XX = np.matmul(X.T,X)/n  # Gram matrix

    Xy = np.matmul(X.T,Y)

    error = Y-np.matmul(X,beta)  # residuals

    Qhat = (error**2).sum()/n  # average of squared residuals



    # Main Loop
    Iter = 0
    while Iter< MaxIter:
        Iter+=1
        beta_old = beta

        ## Go over each coordinate

        for j in range(p):

            # Compute the shoot and update the variable
            S0 = np.matmul(XX[j,:],beta) - (XX[j,j]*beta[j])-Xy[j]
            # \
            if abs(beta[j])>0:

                error = error + (X[:,j]*beta[j])
                Qhat = (error**2).sum()/n


            if n**2<((lam*gamma[j])**2 /XX[j,j]):
                beta[j]=0

            elif S0> ((lam/n)*gamma[j]*np.sqrt(Qhat)):
                beta[j]=(  ( lam* gamma[j] / np.sqrt( n**2 - (lam * gamma[j])**2 / XX[j,j]  ) )  *  np.sqrt( max(Qhat - (S0**2/XX[j,j]),0) )  - S0 )   /   XX[j,j]
                error=error-X[:,j]*beta[j]

            elif S0 < (-1* (lam/n)*gamma[j]*np.sqrt(Qhat)):
                # Optimal beta(j) > 0
                beta[j] = (  -  ( lam * gamma[j] / np.sqrt( n**2 - (lam * gamma[j])**2 / XX[j,j]  ) )  *  np.sqrt( max( Qhat - (S0**2/XX[j,j]),0) )  - S0 )   /   XX[j,j]
                error = error - X[:,j]*beta[j]
            elif abs(S0)<= (lam/n)*gamma[j]*np.sqrt(Qhat):
                beta[j]=0


        #Update the primal and dual

        fobj = np.sqrt( ((np.matmul(X,beta)-Y)**2).sum()/n ) + (np.transpose(lam*gamma/n))* abs(beta)

        if np.linalg.norm(error)>1.0e-10:
            aaa = np.sqrt(n)*error /np.linalg.norm(error)
            dual = aaa.T @ Y/n -abs(lam*gamma/n-abs(np.matmul(X.T,aaa)/n)).T@abs(beta)

        else:
            dual = lam * (gamma.T)@abs(beta)/n


    betaSQ = beta

    return betaSQ
