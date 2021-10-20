import matplotlib.pyplot as plt
import seaborn as sns

def plot_HMscatter(df_HM_Detection, fp_output):
    for i in range(len(df_HM_Detection)):
    df_HM_Detection.HMscore[i] = round(float(df_HM_Detection.HMscore[i]), 3)
    df_HM_Detection.mutNumber[i] = int(df_HM_Detection.mutNumber[i])
    sns.set_style('ticks')
    sns.set_context('paper')
    sns.scatterplot(data=df_HM_Detection, x='mutNumber',
                    y='HMscore', hue='PRmark', style='HMmark')
    plt.axhline(y=cutoff_HMscore, color='black', linestyle='--')
    plt.axvline(x=cutoff_mutLoad, color='black', linestyle='--')
    plt.xscale('log')
    plt.xlabel('Number of mutation')
    plt.ylabel('HM score')
    plt.legend()

    return(0)

def plot_HMFrac():
    return(0)

def plot_HM_SMratio():
    return(0)
