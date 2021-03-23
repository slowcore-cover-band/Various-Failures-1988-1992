##Sign Restrictions
##2/21/20

###Meant to work in tandem with Chris Sims' VAR Tools functions

draw_signs = function(varout, signm, printevery = 50, printfails=5000, ndraws=1000, horiz = 40, bands = c(0.05,0.16,0.5,0.84,0.95), verbose=TRUE){
###Signm has the following structure nv*nshock. Each entry defines the length of the restriction
###and the sign. For example, -5 in the 1,2 position says that the second shock has a negative
###restriction on variable one for five periods.
B = varout$By
nvar = dim(B)[2]
irf_draws = array(0,c(nvar,nvar,horiz,ndraws))
##Setting up S Matrix
bigS = matrix(0,sum(abs(signm)),nvar*max(abs(signm)))
todos = which(signm !=0)
todosi = which(signm !=0,arr.ind=TRUE)
#print(todosi)
#    for (i in 1:length(todos)){
#        if(i>1){
#            prevs = sum(abs(signm[todos[1:i-1]]))
#        }else{
#            prevs = 0
#        }
#        hormaxcur = abs(signm[todos[i]])
#        for(j in 1:hormaxcur){
#            if(signm[todos[i]]>0){
#                bigS[(prevs+j),(nvar*(j-1)+todosi[i,1])] = 1
#            }else{
#                bigS[(prevs+j),(nvar*(j-1)+todosi[i,1])] = -1
#            }
#        }
#    }
    i = 1
    totalcount = 1
    while(i<ndraws){
        ##Drawing VAR parameters
        varpar = postdraw(varout,1)

        ##Drawing Rotation matrix
        rotate = matrix(rnorm(nvar*nvar),nvar,nvar)
        rotate = qr.Q(qr(rotate))

        ##getting new A_0
        Ai = rotate %*% varpar$smat[,,1]

        ##impulse draw
        impulse = impulsdtrf(vout=list(By = varpar$By[,,,1]), smat = Ai, nstep = horiz)
        #print(dim(impulse))
        allgood = rep(0,length(todos))
        for(j in 1:length(todos)){
            choriz = abs(signm[todos[j]])
            #print(choriz)
            if(signm[todos[j]]>0){
                check = sum(as.numeric(impulse[todosi[j,1],todosi[j,2],1:choriz]<0))
                if(check == 0){
                    allgood[j] = 1
                }else{
                    allgood[j] = 0
                }
            }else{
                check = sum(as.numeric(impulse[todosi[j,1],todosi[j,2],1:choriz]>0))
                if(check == 0){
                    allgood[j] = 1
                }else{
                    allgood[j] = 0
                }
            }
        }
        sumgood = sum(allgood)
        if(sumgood==length(todos)){
            irf_draws[,,,i]=impulse
            #print("Accepted!")
            i=i+1
            totalcount = totalcount+1
            if(is.null(printevery) == FALSE){
                if(i%%printevery==1){
                    pctacc = 100*(i/totalcount)
                    print(paste0(pctacc,"% of draws accepted or ",i-1," of ",totalcount))
                }
            }
        }else{
            #print("Rejected!")
            i=i
            totalcount = totalcount+1
            if(verbose==TRUE){
                if(is.null(printfails) == FALSE){
                    if(totalcount%%printfails==0){
                        print(paste(totalcount,"matrices drawn so far and",i-1,"found"))
                    }
                }
            }
        }
    }
    ##Computing bands
    irfbands = array(0,c(nvar,nvar,horiz,length(bands)))
    #print(length(bands))
    for(i in 1:nvar){
    	for(j in 1:nvar){
    		for(k in 1:horiz){
                for(l in 1:length(bands)){
                    irfbands[i,j,k,l]=quantile(irf_draws[i,j,k,],bands[l])
                }
    		}
    	}
    }
    return(list(draws = irf_draws,bands = irfbands))
}
