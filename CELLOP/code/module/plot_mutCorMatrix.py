import matplotlib.pyplot as plt
import seaborn as sns

def plot_mutCM(df_mutCorMatrix, fp_output):
    for i in range(len(df_mutCorMatrix)):
        df_mutCorMatrix.dotSize[i] = round(float(df_mutCorMatrix.dotSize[i]), 3)
    sns.set_style('ticks')
    sns.set_context('paper')
    sns.scatterplot(data=df_mutCorMatrix, x='listA',
                y='listB', hue='dotColor', size='dotSize')
    plt.legend('')
    plt.savefig(fp_output+'Figure2_CoMutation.pdf', dpi=300)
    return(0)