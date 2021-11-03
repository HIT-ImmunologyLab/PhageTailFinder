import numpy as np

def viterbi(trainsition_probability, emission_probability, pi, obs_seq):
    trainsition_probability=np.array(trainsition_probability)
    emission_probability=np.array(emission_probability)
    pi=np.array(pi)
    Row = np.array(trainsition_probability).shape[0]
    Col = len(obs_seq)
    F=np.zeros((Row,Col))
    F[:,0]=pi*np.transpose(emission_probability[:,obs_seq[0]])
    for t in range(1,Col):
        list_max=[]
        for n in range(Row):
            list_x=list(np.array(F[:,t-1])*np.transpose(trainsition_probability[:,n]))
            list_p=[]
            for i in list_x:
                list_p.append(i*10000)
            list_max.append(max(list_p)/10000)
        F[:,t]=np.array(list_max)*np.transpose(emission_probability[:,obs_seq[t]])

    ans = []
    for i in range(len(F[0])):
        if F[0][i] > F[1][i]:
            ans.append(1)
        else:
            ans.append(0)
    return ans, F