import numpy as np
import scipy.stats
import pandas

# Prior distribution parameters
U = 3
V = 1000
EPS = 1e-1

def printResults (d, videos, ms):
    ys = []
    yhats = []
    for video in videos:
        idxs = np.nonzero(d.video == video)[0]
        avgLearningGains = np.mean(d.iloc[idxs].posttest - d.iloc[idxs].pretest)
        ys.append(avgLearningGains)
        yhats.append(ms[video])
    print(scipy.stats.pearsonr(yhats, ys))
    print(scipy.stats.spearmanr(yhats, ys))

def MStep (d, workers, ms, s2s, oldRs, oldMus, useMu, sigmaType, useR):
    rs = {}
    mus = {}
    sigma2s = {}
    for worker in workers:
        idxs = np.nonzero(d.worker == worker)[0]
        mu = 0
        r1 = 0
        r2 = 0
        for idx in idxs:
            video = d.iloc[idx].video
            label = d.iloc[idx].helpOthers
            mu += label - oldRs[worker]*ms[video]
            r1 += (label - oldMus[worker]) * ms[video]
            r2 += ms[video]**2 + s2s[video]
        r = r1/r2
        mu /= len(idxs)
        sigma2 = 0
        for idx in idxs:
            video = d.iloc[idx].video
            label = d.iloc[idx].helpOthers
            sigma2 += (label - mu)**2 - 2*(label - mu)*oldRs[worker]*ms[video] + oldRs[worker]**2*(ms[video]**2 + s2s[video])
        if useMu:
            mus[worker] = mu
        else:
            mus[worker] = 0

        if sigmaType == "pretest":
            sigma2s[worker] = 1./(np.mean(d.iloc[idxs].pretest) + EPS)**.5
        elif sigmaType == "sigma":
            sigma2s[worker] = sigma2
        else:
            sigma2s[worker] = 1

        if useR:
            rs[worker] = r
        else:
            rs[worker] = 1
    return rs, mus, sigma2s

def EStep (d, videos, rs, mus, sigma2s):
    ms = {}
    s2s = {}
    for video in videos:
        s2inv = 1/V**2
        m = U/V**2 
        idxs = np.nonzero(d.video == video)[0]
        for idx in idxs:
            worker = d.iloc[idx].worker
            label = d.iloc[idx].helpOthers
            sigma2 = sigma2s[worker] + EPS
            r = rs[worker]
            s2inv += 1/sigma2
            m += (label - mus[worker])/(r*sigma2)
        s2 = s2inv ** -1
        s2s[video] = s2
        m = s2 * m
        ms[video] = m
    return ms, s2s

def EM (d, useMu, sigmaType, useR):
    workers = np.unique(d.worker)
    videos = np.unique(d.video)
    rs = { worker:1. for worker in workers }
    mus = { worker:0. for worker in workers }
    sigma2s = { worker:1. for worker in workers }
    for i in range(50):
        ms, s2s = EStep(d, videos, rs, mus, sigma2s)
        rs, mus, sigma2s = MStep(d, workers, ms, s2s, rs, mus, useMu, sigmaType, useR)
        printResults(d, videos, ms)
    printResults(d, videos, ms)
    return sigma2s

if __name__ == "__main__":
    d = pandas.read_csv("mturk_data.csv")
    for useMu in [ True, False ]:
        for sigmaType in [ "pretest", "sigma", "" ]:
            for useR in [ True, False ]:
                print("{} {} {}".format(useMu, sigmaType, useR))
                sigma2s = EM(d, useMu, sigmaType, useR)
