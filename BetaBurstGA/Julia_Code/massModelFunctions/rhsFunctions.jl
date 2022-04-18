#################################
# Functions for use in the SDE  #
# solver                        #
#################################

function f(du,u,p,t)
    σE,σI,τxE,τxI,
    κSEE,κSIE,κSEI,κSII,
    αEE,αIE,αEI,αII,
    κVEE,κVIE,κVEI,κVII,
    VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E,η_0I,τE,τI = p

    rE=u[1]
    rI=u[2]
    vE=u[3]
    vI=u[4]
    pEE=u[5]
    gEE=u[6]
    pIE=u[7]
    gIE=u[8]
    pEI=u[9]
    gEI=u[10]
    pII=u[11]
    gII=u[12]
    #rE
    du[1] =(1.0/τE)*(-gEE*rE -gEI*rE - κVEE*rE - κVEI*rE +2.0 * rE * vE + (ΔE / (τE*π)))
    #rI
    du[2] =(1.0/τI)*(-gII*rI - gIE*rI -κVIE*rI - κVII*rI + 2.0 * rI * vI + (ΔI / (τI*π)))
    #vE
    du[3] =(1.0/τE)*(gEE*(VsynEE - vE) + gEI*(VsynEI - vE) + κVEI*(vI - vE) - (τE^2.0)*(π^2.0) * (rE^2) +  vE^2.0 + η_0E + u[13])
    #vI
    du[4] =(1.0/τI)*(gIE*(VsynIE -vI) + gII*(VsynII - vI) + κVIE*(vE - vI) - (τI^2.0)*(π^2.0)* (rI^2) + vI^2.0 + η_0I + u[14])
    #psiEE
    du[5] = αEE * (-pEE + κSEE * rE)
    #gEE
    du[6] = αEE * (-gEE + pEE)
    #pIE
    du[7] = αIE * (-pIE + κSIE * rE)
    #gIE
    du[8] = αIE * (-gIE + pIE)
    #pEI
    du[9] = αEI * (-pEI + κSEI * rI)
    #gEI
    du[10] = αEI * (-gEI + pEI)
    #pII
    du[11] = αII * (-pII + κSII * rI)
    #gII
    du[12] = αII * (-gII + pII)

    du[13] = -u[13]/τxE

    du[14] = -u[14]/τxI


end

function g(du,u,p,t)
    σE,σI,τxE,τxI,
    κSEE,κSIE,κSEI,κSII,
    αEE,αIE,αEI,αII,
    κVEE,κVIE,κVEI,κVII,
    VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E,η_0I,τE,τI = p

    du[1] = 0.
    du[2] = 0.
    du[3] = 0.
    du[4] = 0.
    du[5] = 0.
    du[6] = 0.
    du[7] = 0.
    du[8] = 0.
    du[9] = 0.
    du[10] = 0.
    du[11] = 0.
    du[12] = 0.
    du[13] = σE
    du[14] = σI
    
end

function fdeterministic(du,u,p,t)
    κSEE,κSIE,κSEI,κSII,
    αEE,αIE,αEI,αII,
    κVEE,κVIE,κVEI,κVII,
    VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E,η_0I,τE,τI = p

    rE=u[1]
    rI=u[2]
    vE=u[3]
    vI=u[4]
    pEE=u[5]
    gEE=u[6]
    pIE=u[7]
    gIE=u[8]
    pEI=u[9]
    gEI=u[10]
    pII=u[11]
    gII=u[12]
    #rE
    du[1] =(1.0/τE)*(-gEE*rE -gEI*rE - κVEE*rE - κVEI*rE +2.0 * rE * vE + (ΔE / (τE*pi)))
    #rI
    du[2] =(1.0/τI)*(-gII*rI - gIE*rI -κVIE*rI - κVII*rI + 2.0 * rI * vI + (ΔI / (τI*pi)))
    #vE
    du[3] =(1.0/τE)*(gEE*(VsynEE - vE) + gEI*(VsynEI - vE) + κVEI*(vI - vE) - (τE^2.0)*(pi^2.) * (rE^2.0) +  vE^2.0 + η_0E )
    #vI
    du[4] =(1.0/τI)*(gIE*(VsynIE -vI) + gII*(VsynII - vI) + κVIE*(vE - vI) - (τI^2.0)*(pi^2.)* (rI^2.0) + vI^2.0 + η_0I )
    #pEE
    du[5] = αEE * (-pEE + κSEE * rE)
    #gEE
    du[6] = αEE * (-gEE + pEE)
    #pIE
    du[7] = αIE * (-pIE + κSIE * rE)
    #gIE
    du[8] = αIE * (-gIE + pIE)
    #pEI
    du[9] = αEI * (-pEI + κSEI * rI)
    #gEI
    du[10] = αEI * (-gEI + pEI)
    #pII
    du[11] = αII * (-pII + κSII * rI)
    #gII
    du[12] = αII * (-gII + pII)

end

function ftest(rE,rI,vE,vI,pEE,gEE,pIE,gIE,pEI,gEI,pII,gII,p)
    κSEE,κSIE,κSEI,κSII,
    αEE,αIE,αEI,αII,
    κVEE,κVIE,κVEI,κVII,
    VsynEE,VsynIE,VsynEI,VsynII,ΔE,ΔI,η_0E,η_0I,τE,τI = p

   
    #rE
    return [(1/τE)*(-gEE*rE -gEI*rE - κVEE*rE - κVEI*rE +2 * rE * vE + (ΔE / (τE*pi))),
    
    (1/τI)*(-gII*rI - gIE*rI -κVIE*rI - κVII*rI + 2 * rI * vI + (ΔI / (τI*pi))),
    
    (1/τE)*(gEE*(VsynEE - vE) + gEI*(VsynEI - vE) + κVEI*(vI - vE) - (τE^2)*(pi^2) * (rE^2) +  vE^2 + η_0E ),
    
    (1/τI)*(gIE*(VsynIE -vI) + gII*(VsynII - vI) + κVIE*(vE - vI) - (τI^2)*(pi^2)* (rI^2) + vI^2 + η_0I ),
    
    αEE * (-pEE + κSEE * rE),
    
    αEE * (-gEE + pEE),
    
    αIE * (-pIE + κSIE * rE),
    
    αIE * (-gIE + pIE),
    
    αEI * (-pEI + κSEI * rI),
    
    αEI * (-gEI + pEI),
    
    αII * (-pII + κSII * rI),
   
     αII * (-gII + pII)]

end
