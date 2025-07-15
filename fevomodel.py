## ---------------------------
##
## Script name: EcoEvoModel.py
##
## Purpose of script: This code runs the model and generates simulations.
## the model is an eco-evolutionary dynamics of community assembly, where
## new species are added to the system as the dynamics run. The ecology
## is determined by a differential equation of species interactions.
## If you use our generated dataframes, this script is not necessary.
##
## Author: gui araujo
## Computational Ecology Lab - Department of Biosciences
## Swansea University, UK
## 
##
## Date Created: 2024
##
## Email: gui.dav.araujo@gmail.com
##
## ---------------------------


# Packages used
import numpy as np
from scipy.integrate import odeint
import pickle

# this is the path to the file's directory
import os 
dir_path = os.path.dirname(os.path.realpath(__file__))


def IncludeInteraction(c0,c1,a,b):
# Goal: create an interaction by generating weights in the interaction matrix
# Inputs: indexes of species to receive interaction (c0,c1), matrix of positive interactions (a), matrix of negative interactions (b)
# Outputs: matrix of positive interactions (a) and matrix of negative interactions (b), modified
        
    # strengths are drawn
    w1 = sigma*np.random.randn()
    w2 = sigma*np.random.randn()
    
    
    # for propInvasion model, the type of interaction is chosen according to the drawn proportions
    if(allowInvasions==1 and propInvasions==1):
        
        piChoice = np.random.rand()
        if(piChoice<piM):
            piChosen='M'
        elif(piChoice>piM and piChoice<piM+piP):
            piChosen='P'
        else:
            piChosen='C'
            
        if(piChosen=='M'):
            w1 = abs(sigma*np.random.randn())
            w2 = abs(sigma*np.random.randn())
        elif(piChosen=='C'):
            w1 = -abs(sigma*np.random.randn())
            w2 = -abs(sigma*np.random.randn())
        elif(piChosen=='P'):
            if(np.random.rand()<0.5):
                w1 = abs(sigma*np.random.randn())
                w2 = -abs(sigma*np.random.randn())
            else:
                w1 = -abs(sigma*np.random.randn())
                w2 = abs(sigma*np.random.randn())

    # encoding strenghts in the positive and negative matrices
    if(w1>0):
        a[c0][c1]=abs(w1)
    else:
        b[c0][c1]=abs(w1)
    if(w2>0):
        a[c1][c0]=abs(w2)
    else:
        b[c1][c0]=abs(w2)
    


            
            
    # if interaction is of consumer-resource type, ensure that the benefit of the consumer is not larger than the decrease of resource
    if(a[c0][c1]!=0 and b[c1][c0]!=0):
        if(a[c0][c1]>b[c1][c0]):
            a[c0][c1]=b[c1][c0]        

    if(a[c1][c0]!=0 and b[c0][c1]!=0):
        if(a[c1][c0]>b[c0][c1]):
            a[c1][c0]=b[c0][c1]    
            


            
    return a,b





def System(yy,t,pos,neg,rl,dl,ins,S, PP, MM, pp, CC,cl):
# Goal: calculate the differentials of the ecological equation
# Inputs: vector of abundances, time (not used in the equation), matrix of positive interactions (not being used),
# matrix of negative interactions (not being used), vector of growth rates, matrix of intraspecific competition,
# total size of the community (not being used), matrix of mutualisms, matrix of consumers, matrix of resources, matrix of competition, vector of total harvesting costs
# Outputs: vector of differentials    

    # denominators of type II functional response
    p_load = 1/(1 + hp*np.dot(PP,yy))      
    mp_load = 1/(1+hm*np.dot(MM, yy))
             
    # ecological equation
    dy = yy*np.dot(PP,yy)*p_load - yy*np.dot(pp,np.multiply(yy,p_load)) - yy*np.dot(CC,yy) + yy*np.dot(MM,yy)*mp_load +(rl-cl)*yy- yy*np.dot(ins,yy)
        

    
    return dy





def Dynamics(y0,pos,neg,rl,dl,ins,exth,time):
# Goal: run the entire simulation, including all the assembly, outputting the final community
# Inputs: initial vector of abundances, initialised matrices and vectors, extinction threshold, maximum simulation time
# Outputs: several quantities, but most importantly the histories of densities, interaction matrices, growth rates, intraspecific competitions   
 
    # initialise several lists and quantities to help and log the process
    TT = time
    densities = []
    stepcounter = 0
    finish=False
    proceed=False
    SC = S0
    prev=np.copy(y0)
    extinct=[]
    n_ext=0
    n_evo=0
    n_attempts=0
    n_attemptsFull=0
    m_snaps=[]
    ins_snaps=[]
    y_snaps=[]
    new_snaps=[]
    rk_list=[]
    target_snaps=[]
    attempt_snaps=[]
    new=0
    target=0
    movingDensities=[]
    isOscillating=0
    
    overrideEq=0
    isInvasion=0
    
    n_created=0
    n_destroyed=0
    

    # generate matrices for each interaction types, using pos and neg (in the final version, these all start out as zeros)
    aux = np.zeros_like(pos)
    aux[pos != 0] = 1
    aux = np.logical_and(aux, aux.T)
    PP = np.where(aux == 0, pos, 0)
    MM = np.where(aux != 0, pos, 0)
    
    bux = np.zeros_like(neg)
    bux[neg != 0] = 1
    bux = np.logical_and(bux, bux.T)
    pp = np.where(bux == 0, neg, 0)
    CC = np.where(bux != 0, neg, 0)
    
    # vector of harvesting costs (also zero at the beginning)
    cl = cost * ( np.count_nonzero(PP,axis=1)+np.count_nonzero(MM,axis=1) )


    # loops through the entire size of the simulation time
    for iii in range(TT):
        
        # calls the numerical integration using the equation in System
        # it integrates only until it reaches a probetime
        step = odeint(System, y0.copy(), np.linspace(0,probetime,int(probetime/timestep)), args=(pos,neg,rl,dl,ins,SC,PP,MM,pp,CC,cl))

        # collects abundances
        y0 = step[-1].copy()
        
        
        # this block is to check who's extinct and store alive indexes
        # isolated species are removed after 20 assembly events (they stay in the beginning to give time to establish interactions)
        alive=[]
        alive_weights=[]
        for k in range(SC):

            if(k not in extinct):

                if(y0[k]<exth or (sum(pos[k,:]+pos[:,k]+neg[k,:]+neg[:,k])==0 and n_evo>20)):
                    
                    if(k not in extinct):
                        
                        if(sum(pos[k,:]+pos[:,k]+neg[k,:]+neg[:,k])==0):
                            print('!!!DISCONNECTED')
   
                        y0[k]=0
                        pos[k, :] = 0
                        pos[:, k] = 0
                        neg[k, :] = 0
                        neg[:, k] = 0
                        ins[k, :] = 0
                        ins[:, k] = 0
                        rl[k]=0
                        dl[k]=0

                        extinct.append(k)
                        n_ext+=1
                        print('!!!EXTINCT: '+str(k))
                        
                else:
                    # there's a limit of abundances. If reached, then the simulation blows up and stops
                    if(y0[k]>100000000000):
                        print('DIVERGENCE!')
                        stop
                    if(k not in alive):
                        alive.append(k)
                        alive_weights.append(y0[k])



        # this counts the attempts to assembly, given proposals of new species
        evo_attempts=0
        


        if(len(alive)>1):
            
            evolve=True


            # this code checks for equilibrium to trigger an assembly event.
            if(overrideEq<override_size and (overrideEq<=2 or isOscillating==0)): # if the system is oscillating, then it keeps going until overriden by an assembly
                for k in range(SC):
                    
                    if(y0[k]>0):
                        if(np.abs(prev[k]-y0[k])/y0[k]>0.0001): # equilibrium is defined as relative changes in all abundances being less than 0.01% between probes
                            evolve=False


            
            
            prev=np.copy(y0) # store current densities to serve as previous in the next iteration, to check the equilibrium
            
            
            
            if(evolve): # if an assembly event occurs, this block runs
                
                
                
                if(n_evo>=n_evolutions): # this marks the final assembly event to trigger the end of the simulation
                    finish=True
                    print(str(SC)+','+'))))))))))))))))))))))))): '+str(n_evo)+' evolved')
                else:
                
                    # if there's an index available, the new species occupies it
                    if(len(extinct)>0):
                        new = extinct[0]
                        del extinct[0]
                    else:
                        new=SC
                        SC+=1 # if not, the network will expand
                        
                        # this block is to restart all arrays expanded to accomodate a new index
                        rs = np.zeros(SC)
                        ds = np.zeros(SC)
                        ys = np.zeros(SC)
                        inss = np.zeros((SC,SC))
                        poss = np.zeros((SC,SC))
                        negs = np.zeros((SC,SC))
                        prevs = np.zeros(SC)
                        
                        for i in range(SC-1):
                            rs[i]=rl[i]
                            ds[i]=dl[i]
                            
                            ys[i]=y0[i]
                            prevs[i]=prev[i]
                            
                            for j in range(SC-1):
                                
                                poss[i][j]=pos[i][j]
                                negs[i][j]=neg[i][j]
                                inss[i][j]=ins[i][j]
                        
                        rl = rs.copy()
                        dl = ds.copy()
                        ins = inss.copy()
                        pos = poss.copy()
                        neg = negs.copy()
                        y0 = ys.copy()
                        prev = prevs.copy()
                            
                    
                    
                    
                            
                    # this loop will repeat until a proposed species is accepted (because it can actually grow)
                    while(proceed==False):
                        
                        if(allowInvasions):
                            isInvasion=1
                        
                        
                        # this will determine how incomplete the inheritance will be, based on nicheshift (Delta), to be used by evolution
                        shift = np.random.randint(1,nicheshift+1)
                        n_created = np.random.randint(shift+1)
                        n_destroyed = shift - n_created


                                
                                
                        evo_attempts+=1

                        # this chooses a parent species, amongst the living ones
                        chosen=np.random.choice(alive)

                            
                        
                        n_attempts+=1
                        # call the function to include a new species
                        y0,pos,neg,rl,dl,ins,target = EvoSpecies(y0,pos,neg,rl,dl,ins,chosen,SC,extinct,new,isInvasion,n_created,n_destroyed)
                        
                        
                        # recalculate matrices for each interaction type and the cost
                        aux = np.zeros_like(pos)
                        aux[pos != 0] = 1
                        aux = np.logical_and(aux, aux.T)
                        PP = np.where(aux == 0, pos, 0)
                        MM = np.where(aux != 0, pos, 0)
                        
                        bux = np.zeros_like(neg)
                        bux[neg != 0] = 1
                        bux = np.logical_and(bux, bux.T)
                        pp = np.where(bux == 0, neg, 0)
                        CC = np.where(bux != 0, neg, 0)
                        
                        cl = cost * ( np.count_nonzero(PP,axis=1)+np.count_nonzero(MM,axis=1) )


                        # calculate the differentials to see if the new species will grow
                        dy = System(y0,0,pos,neg,rl,dl,ins,SC,PP,MM,pp,CC,cl)
                        
                        if(dy[new]>0):# if it grows, then proceed
                            
                            proceed=True


                            n_evo+=1
                            
                            
                            # since the assembly was valid, this will log a bunch of values
                            if(evo_attempts>1):
                                print('>>>>>>>>>>>>>>>>>> ATTEMPTS: '+str(evo_attempts))
                                attempt_snaps.append(np.array([isInvasion,evo_attempts]))
                            if(isInvasion==0):
                                print('-------------- Tgt: '+str(target))
                            if(isInvasion==1):
                                print('$$$$$$$$$$ INVASION $$$$$$$$$')
                            
                            print('new:'+str(new)+' sc:'+ str(SC) +' ext:'+str(len(extinct)))
                            print('MPC: '+str(np.count_nonzero(MM)/2)+', '+str(np.count_nonzero(PP))+', '+str(np.count_nonzero(CC)/2) +' - dens: '+ str(np.round(np.sum(y0),2)))
                            print('>>>>>NICHESHIFT:',shift)
                            print('***** EVO: '+str(n_evo))
                            
                        # if the new species doesn't grow, all the changes are undone
                        else:
                            if(evo_attempts<max_attempts):
                                y0[new]=0
                                for i in range(SC):
                                    pos[new][i] = 0
                                    pos[i][new] = 0
                                    neg[new][i] = 0
                                    neg[i][new] = 0
                                    ins[new][i]=0
                                    ins[i][new]=0
                                rl[new]=0
                                dl[new]=0

                        
                            if(evo_attempts>=max_attempts):
                                y0[new]=0
                                for i in range(SC):
                                    pos[new][i] = 0
                                    pos[i][new] = 0
                                    neg[new][i] = 0
                                    neg[i][new] = 0
                                    ins[new][i]=0
                                    ins[i][new]=0
                                rl[new]=0
                                dl[new]=0
                                
                                
                                if(isInvasion==1):
                                    print('$$$$$ INVASION $$$$$')
                                print('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*', max_attempts,' ATTEMPTS')
                                print('-------------- Tgt: '+str(target))
                                print('new:'+str(new)+' sc:'+ str(SC) +' ext:'+str(len(extinct)))
                                print('MPC: '+str(np.count_nonzero(MM)/2)+', '+str(np.count_nonzero(PP))+', '+str(np.count_nonzero(CC)/2) +' - dens: '+ str(np.round(np.sum(y0),2)))

                                print('******** EVO: '+str(n_evo))
                                print('*-*-*-*-*-*-*-*-*-*-*-*-*-*-*')

                                proceed=True
                                if(n_evo>50):
                                    finish=True # after 50 assembly events, a total of max_attempts (5000) results in the termination of the simulation [this is not supposed to happen, it's just a lock]

                                
            
        if(proceed):
            soly = np.zeros(SC) # the official vector of abundances is updated
            for k in range(SC):
                soly[k]=y0[k]
            
            

            # in this block, variables are prepared to be saved
            # the matrices are transformed into a compact form that holds indexes and nonzero values
            mat = pos-neg
            lin = [np.array([mat.shape[0],0,0])]

            for ii in range(mat.shape[0]):
                for jj in range(mat.shape[1]):
                    if(mat[ii][jj]!=0):
                        lin.append(np.array([ii,jj,mat[ii][jj]]))

            m_snaps.append(np.array(lin))

            y_snaps.append(np.array([soly]).T)
            new_snaps.append(new)
            target_snaps.append(target)
            
            proceed=False
            overrideEq=0
            isInvasion=0
            if(len(movingDensities)==0):
                isOscillating=0
            movingDensities=[]
            
        
        
        stepcounter+=1
        
        # this block probes if there are oscillations occuring
        noise=0
        isOscillating=0
        if(stepcounter%100==0 and overrideEq!=0):
            movingDensities.append(np.array([y0]).T)
            noise = np.divide(np.array(movingDensities).std(axis=0), np.array([y0]).T, out=np.zeros_like(np.array([y0]).T), where=np.array([y0]).T!=0).sum()/len(alive)
            if(noise>1):
                isOscillating=1
            else:
                isOscillating=0
        
        if(stepcounter%5000==0):
            if(len(movingDensities)!=0):
                print(stepcounter/100, np.round(noise,3), 'OSC=',isOscillating)
            else:
                print(stepcounter/100, 0, 'OSC=',isOscillating)
            overrideEq+=1

            
        if(len(alive)<=1):
            finish=True # if there's less than 2 species alive, the simulation is terminated
            
                

        
        if(finish): # this breaks out of the simulation loop
            break
       
    # this is outside of the loop, after the simulation ended
    soly = np.zeros(SC)
    for k in range(SC):
        soly[k]=y0[k]
       
    densities.append(soly)


    rk_list.append([rl,ins,soly,pos,neg])


    return np.array(densities), [n_ext,n_evo,np.array(m_snaps, dtype=object),np.array(y_snaps, dtype=object),n_attempts,n_attemptsFull,0,0,np.array(new_snaps, dtype=object),np.array(target_snaps, dtype=object),np.array(attempt_snaps, dtype=object),np.array(ins_snaps, dtype=object),0,isOscillating], rk_list




def EvoSpecies(y0,a,b,rl,dl,ins,target,SC,extinct,new,isInvasion,n_created,n_destroyed):
# Goal: propose a new species for an assembly event
# Inputs: vector of abundances y0, matrix of positive interactions a, matrix of negative interactions b,
# independent growth rates rl, death rates dl, matrix of intraspecific competition ins,
# parent species index target (for evo), network size SC, list of extinct entries extinct, new species index new,
# bool for invasion isInvasion, number of created interactions n_created (for evo), number of destrpyed interactions n_destroyed (for evo)
# Outputs: new values for vector of abundances y0, matrix of positive interactions a, matrix of negative interactions b,
# independent growth rates rl, death rates dl, matrix of intraspecific competition ins, and the parent species target (no purpose anymore)
    
    global piM
    global piP
    global piC
    
    y0[new]=exth # new species is initialised at the extinction threshold
   
    
    # initialise new variables
    rl[new] = mur+ mur*sigr*(np.random.randn())
    dl[new] = 0
    ins[new][new] = rl[new]/np.random.lognormal(mui,sigi)
   
    # invasion is simple, just randomly choose interactions
    if(isInvasion==1):
        
        if(propInvasions==1): # Choosing proportions if propInvasion model is selected
            piM=np.random.rand()
            piP=np.random.rand()
            piC=np.random.rand()
            piT=piM+piP+piC
            piM=piM/piT
            piP=piP/piT
            piC=piC/piT
        
        # choose the connectance of the new species based on p1 and p2
        C = p1+(p2-p1)*np.random.rand()
    
        # include the interactions in a and b
        for i in range(SC):
            if(i!=new and y0[i]>0):
                
                if(np.random.rand()<C):
                    a,b = IncludeInteraction(new, i, a, b)
                
                        
        
    # this is for evolution
    if(isInvasion==0):
        
        
    
        nonzero=[] # list to store the indexes of species with interactions with the parent
        zero=[] # list to store indexes wihtout interactions with the parent
        
        for i in range(SC):
            if(i!=new and i!=target and y0[i]>0):
    
                # in this block, all interactions of the parent are inherited
                if(b[target][i]>0):

                    b[new][i] = b[target][i]*np.abs(1 +msig*np.random.randn())
                    
                        
                if(b[i][target]>0):

                    b[i][new] = b[i][target]*np.abs(1 +msig*np.random.randn())
                    
                    
                    
                if(a[target][i]>0):
                    a[new][i] = a[target][i]*np.abs(1 +msig*np.random.randn())


                    if(b[i][new]>0): # this means that new is a consumer of i
                        if(a[new][i]>b[i][new]): # then it can't be larger
                            a[new][i]=b[i][new]
                    
    
                if(a[i][target]>0):

                    a[i][new] = a[i][target]*np.abs(1 +msig*np.random.randn())
                    
                    if(b[new][i]>0): # this means that i is a consumer of new
                        if(a[i][new]>b[new][i]): # then it can't be larger
                            a[i][new]=b[new][i]
                        

                    
                    
                
                # if there is any interaction, i goes to nonzero. If not, i goes to zero
                if(a[i][new]>0 or a[new][i]>0 or b[i][new]>0 or b[new][i]>0):    
                    if i not in nonzero:
                        nonzero.append(i)
                else:
                    zero.append(i)

        
        # the number of new interactions created and old interactions destroyed is already defined (using nicheshift, which is lambda)
        toCreate = n_created
        toDestroy = n_destroyed
        
        # all indexes in nonzero are eligible to destruction, all indexes in zero are eligible to creation
        # if there are not enough interactions to destroy, all but one will be destroyed
        # if there are no available non-interactions for creation, all will be created
        if(len(nonzero)!=0):
            if(toDestroy>len(nonzero)):
                if(len(nonzero)>1):
                    destroy = nonzero
                else:
                    destroy=[]
            else:
                destroy = np.random.choice(nonzero,size=toDestroy,replace=False)
            
            for i in destroy:
                a[i][new]=0
                a[new][i]=0
                b[i][new]=0
                b[new][i]=0
        
        
        if(len(zero)!=0):
            if(toCreate>len(zero)):
                n_create = len(zero)
            else:
                n_create = toCreate
                
            
            for i in range(n_create):


                i = np.random.choice(zero)
                zero.remove(i)

                a,b = IncludeInteraction(new, i, a, b)
                

    return y0,a,b,rl,dl,ins,target

    
  
####################################################### 
      


timestep=0.01 # scale of timesteps for numerical integration
probetime=10*timestep # interval of integration before checking for equilibrium to include an assembly event
override_size=1 # how long to wait before continuing the assembly when oscillations are detected
max_attempts=5000 # maximum number of trials for sampling of a thriving new species (never reached)



# final time of simulation (not reached, because it stops based on assembly events)
ttime = 500000*100


S0=5 # Initial number of species


# parameters of normal for independent growth rate
mur = 0.1
sigr = 0.1

# parameters of lognormal for intraspecific competition
mui = -2.2
sigi = 0.5

msig = 0.05 # std of the change in weights when evolving
exth = 0.000001 # extinction threshold

hm=0.1 # mutualism handling time for type II
hp=0.1 # consumer-resource handling time for type II


########################## PARAMETERS TO CHANGE

allowInvasions=0 # This controls whether it's evolution or invasions. Set 0 for evolution and 1 for invasions
p1=0.01
p2=0.5
propInvasions=1 # Change to control the type of invasion. 0 for random interactions and 1 for random proportions of interaction types for the invader species
# global variables for controlling proportions if propInvasions=1 (not to change)
piM=0
piP=0
piC=0

nicheshift=5 # This is the Delta (how incomplete is inheritance). Change it to 5 or 30 to reproduce strong and weak evolution

sigma=0.05 # This is the sigma of interactions. Change it between [0.01,0.2] to reproduce results


costLevel=5 # the costs we used range from level 1 to 6
cost=costLevel*0.5*0.01 # This is the lambda (harvesting cost). Change it between [0.005,0.03] to reproduce results, but better to change by levels instead


num_sample=1 # code for the sample, change it to save different numbers of samples
print('{{{{{{{{{{{{{{{{{ num_sample= '+str(num_sample))

n_evolutions = 500 # Number of assembly events

#######################


# RUN
isTest=0 # if running a new simulation with new parameter values chosen above, not exactly the batches used in the paper, use isTest=1. Otherwise, use isTest=0

isHeatmap=0 # if isTest=0: isHeatmap=1 runs all combinations of simulations used in the heatmaps. Otherwise, if isHeatmap=0, it runs all 8 scenarios for other figures (don't change any parameters for these)




# This will run custom simulations
if(isTest==1):
    print('{{{{{{{{{{{{{{{{{ num_sample= '+str(num_sample))

    sols = []
    
    
    a = np.zeros((S0,S0)) # matrix of positive interactins
    b = np.zeros((S0,S0)) # matrix of negative interactins
    r = np.zeros(S0) # vector of independent growth rates
    d = np.zeros(S0) # vector of death rates [not being used]
    ins = np.zeros((S0,S0)) # matrix of intraspecific competition
    y0 = np.zeros(S0) # vector of abundances
     
    for i in range(S0):
           
        r[i] = mur+ mur*sigr*(np.random.randn())
        d[i] = 0
        ins[i][i] = r[i]/np.random.lognormal(mui,sigi)
    
    
    for k in range(S0):
        y0[k] = 0.02*np.random.rand()
    
            
    # run and store simulation
    densities, numbers, rk_list = Dynamics(y0,a,b,r,d,ins,exth,ttime)
    sols.append((densities, numbers, rk_list))

    # save results
    #path = dir_path+'/data/modelData_1/test_1/'+str(num_sample)+'FEVO-O1.pkl'
    
    #with open(path, 'wb') as file:
          
    #    pickle.dump(sols, file)
        
        
        
    # save results
    path = dir_path+'/data/analysisData_1/'+str(num_sample)+'FEVO-O1.pkl' # Change this path to where you want to save the simulation
    
    with open(path, 'wb') as file:
          
        pickle.dump(sols, file)




# This will run the 8 scenarios for figures 2, 3, and 4
if (isTest==0 and isHeatmap==0):
    
    for community_type in ['t1','t2']: # runs for both community types
        
        if(community_type=='t1'): # defines T1
            sigma=0.05
            cost=0.025   
        if(community_type=='t2'): # defines T2
            sigma=0.2
            cost=0.01
    
        for model in ['evo5','evo30','invR','invP']: # each of the 4 scenarios
        
        
            if(model=='evo5'): # define Evo. Delta=5
                allowInvasions=0
                nicheshift=5
            if(model=='evo30'): # define Evo. Delta=30
                allowInvasions=0
                nicheshift=30
            if(model=='invR'): # define Inv. Rand.
                allowInvasions=1
                propInvasions=0
            if(model=='invP'): # define Inv. Prop.
                allowInvasions=1
                propInvasions=1
                
            
            for sample in range(1,21): # runs and saves 20 samples for each scenario
                num_sample=sample
                print('{{{{{{{{{{{{{{{{{ num_sample= '+str(num_sample))
            
                sols = []
                
                
                a = np.zeros((S0,S0)) # matrix of positive interactins
                b = np.zeros((S0,S0)) # matrix of negative interactins
                r = np.zeros(S0) # vector of independent growth rates
                d = np.zeros(S0) # vector of death rates [not being used]
                ins = np.zeros((S0,S0)) # matrix of intraspecific competition
                y0 = np.zeros(S0) # vector of abundances
                 
                for i in range(S0):
                       
                    r[i] = mur+ mur*sigr*(np.random.randn())
                    d[i] = 0
                    ins[i][i] = r[i]/np.random.lognormal(mui,sigi)
                
                
                for k in range(S0):
                    y0[k] = 0.02*np.random.rand()
                
                        
                # run and store simulation
                densities, numbers, rk_list = Dynamics(y0,a,b,r,d,ins,exth,ttime)
                sols.append((densities, numbers, rk_list))
                
                
                
                # save results
                path = dir_path+'/data/modelData_1/ecoevo_'+community_type+'_'+model+'_1/'+str(num_sample)+'FEVO-O1.pkl'
                
                with open(path, 'wb') as file:
                      
                    pickle.dump(sols, file)
        

# This will run all the combinations of sumulations to create the heatmaps
if (isTest==0 and isHeatmap==1):
    
    # Heatmaps used Evo. Delta=5
    allowInvasions=0
    allowInvasions=1
    nicheshift=5
    
    # This runs and sorts (in a kinda convoluted way) all the samples for each sigma,delta combinations
    for batch in ['','2','3']: # samples were taken in 3 batches of 5 = 15 samples
        for sig in [0.01,0.05,0.1,0.15,0.2,0.25]: # this chooses sigma values
            for clevel in [1,2,3,4,5,6]: # this chooses delta level
                for sample in [1,2,3,4,5]: # this runs 5 samples for a chosen sigma,delta pair
                    
                    sigma=sig
                    cost=clevel*0.5*0.01
                    
                    # this is a bit cumbersome, but it's the way the data was separated into folders before automatisation of the process
                    step=1
                    stepName=''
                    if(clevel>3):
                        stepName='2'
                        step=2
                    
                    sample_multiplier=clevel-1
                    if(step==2):
                        sample_multiplier-=3
                
                    num_sample = int(sample_multiplier)*5+sample
                    
                

                    print('{{{{{{{{{{{{{{{{{ num_sample= '+str(num_sample))
                
                    sols = []
                    
                    
                    a = np.zeros((S0,S0)) # matrix of positive interactions
                    b = np.zeros((S0,S0)) # matrix of negative interactions
                    r = np.zeros(S0) # vector of independent growth rates
                    d = np.zeros(S0) # vector of death rates [not being used]
                    ins = np.zeros((S0,S0)) # matrix of intraspecific competition
                    y0 = np.zeros(S0) # vector of abundances
                     
                    for i in range(S0):
                           
                        r[i] = mur+ mur*sigr*(np.random.randn())
                        d[i] = 0
                        ins[i][i] = r[i]/np.random.lognormal(mui,sigi)
                    
                    
                    for k in range(S0):
                        y0[k] = 0.02*np.random.rand()
                    
                            
                    # run and store simulation
                    densities, numbers, rk_list = Dynamics(y0,a,b,r,d,ins,exth,ttime)
                    sols.append((densities, numbers, rk_list))
                    
                    
                    
                    # save results
                    path = dir_path+'/data/heatData_1/'+batch+'heat'+stepName+'_'+str(sig)+'/'+str(num_sample)+'FEVO-O1'+str(sigma)+'_'+str(clevel)+'.pkl'
                    with open(path, 'wb') as file:
                          
                        pickle.dump(sols, file)
            
            
            
        
