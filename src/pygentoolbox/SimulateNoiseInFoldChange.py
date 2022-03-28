
def main(irslist=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], depth=20, noiselevel=0.5, iterations=100):
    # noiselevel is a float (frequency). Will cause loss of this freq of fragments from total depth of IES.
    # noiselevel is low @ 0 (no noise), high at 1 (complete loss of all fragments due to noise).
    import random
    import seaborn as sns
    import matplotlib.pyplot as plt
    import math

    allnoiseirs, allnoisefc = [], []
    fclist = []
    for irs in irslist:
        ret = int(depth * irs)
        exc = int(depth - ret)
        lfc = math.log2(ret / exc)
        fclist.append(lfc)
        irsnoise, fcnoise = [], []
        for i in range(iterations):
            ies = [0] * exc + [1] * ret  # 1 is retained fragment 0 is excised
            for j in range(int(depth*noiselevel)):
                ies.pop(random.randint(0, len(ies) - 1))
            # collect noisy scores
            newret = sum(ies)
            newexc = len(ies) - newret
            irsnoise.append(newret / len(ies))
            if newexc != 0:
                newfc = newret / newexc
                if newfc != 0:
                    fcnoise.append(math.log2(newfc))
                else:
                    fcnoise.append(lfc)
            else:
                # avoid div by zero
                if newfc != 0:
                    newfc = newret / 1
                else:
                    fcnoise.append(lfc)
        allnoiseirs.append(irsnoise)
        allnoisefc.append(fcnoise)

    names = [str(irs) for irs in irslist]
    namesfc = [str(round(fc, 1)) for fc in fclist]

    fig, axes = plt.subplots(1, 2, figsize=(8, 4))
    sns.boxplot(data=allnoiseirs, ax=axes[0])  # violinplot
    # sns.stripplot(data=allnoiseirs, ax=axes[0], jitter=True)
    axes[0].set(xticklabels=names, xlabel='IRS', ylabel='Noisy IRS')

    sns.boxplot(data=allnoisefc, ax=axes[1])  # violinplot
    # sns.stripplot(data=allnoisefc, ax=axes[1], jitter=True)
    axes[1].set(xticklabels=namesfc, xlabel='Log2FC', ylabel='Noisy Log2FC')
    plt.show()
    plt.close()










