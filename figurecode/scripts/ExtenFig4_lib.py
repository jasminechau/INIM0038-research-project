import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import os



def plot_run_success(success_data):
    success_data = success_data.groupby(['Method', 'Output', 'Sample_groups']).count().reset_index().rename(columns={'Patient': 'Number of Sims'})
    groups = {'low': 'Low', 'med': 'Medium', 'high': 'High'}
    plt.rcParams.update({'font.size': 20, 'xtick.bottom': False, 'ytick.left': False})
    fig, ax = plt.subplots(1,len(groups), sharex=False, sharey=True, figsize=(12*len(groups)+1, 6))
    for j in range(len(groups)):
        use_data = success_data[(success_data.Sample_groups == list(groups.keys())[j])].drop(['Sample_groups'], axis=1)
        plotting_data = use_data.pivot_table(index='Method', columns='Output', values='Number of Sims').fillna(0).reindex(['PhyloWGS', 'CITUP', 'LICHeE', 'Pyclone', 'TRACERx'])
        if 'FAILED' not in plotting_data.columns:
            plotting_data['FAILED'] = 0
        if 'TIMEOUT' not in plotting_data.columns:
            plotting_data['TIMEOUT'] = 0
        if 'COMPLETE' not in plotting_data.columns:
            plotting_data['COMPLETE'] = 0
        plotting_data = plotting_data.fillna(0)
        p1 = ax[j].barh(['PhyloWGS', 'CITUP', 'LICHeE', 'Pyclone', 'TRACERx'], plotting_data['COMPLETE'], 0.9, color=sns.color_palette('deep')[4::-1], linewidth=1, edgecolor='white')
        p2 = ax[j].barh(['PhyloWGS', 'CITUP', 'LICHeE', 'Pyclone', 'TRACERx'], plotting_data['FAILED'], 0.9, left=plotting_data['COMPLETE'], color=sns.color_palette('deep')[4::-1], hatch='x', alpha=0.7, linewidth=1, edgecolor='white')
        p3 = ax[j].barh(['PhyloWGS', 'CITUP', 'LICHeE', 'Pyclone', 'TRACERx'], plotting_data['TIMEOUT'], 0.9, left=np.array(plotting_data['COMPLETE']) + np.array(plotting_data['FAILED']), color=sns.color_palette('deep')[4::-1], hatch='//', alpha=0.3, linewidth=1, edgecolor='white')
        ax[j].set_xlabel('Number of Simulations', fontsize=26)
        ax[j].set_xlim(0, 50)
        ax[j].set_ylabel('', fontsize=26)
        ax[j].set_title('Sample Group: {}'.format(groups[list(groups.keys())[j]]), fontdict={'fontsize':28})
        if j != 0:
            ax[j].set_ylabel('')
        else:
            ax[j].set_yticks(['PhyloWGS', 'CITUP', 'LICHeE', 'Pyclone', 'TRACERx'], labels=['PhyloWGS', 'CITUP', 'LICHeE', 'Pyclone +\nTRACERx', 'TRACERx'])

    complete = mpatches.Patch(facecolor='#636366', label='Complete', edgecolor='white')
    failed = mpatches.Patch(facecolor='#8B898A', label='Failed', edgecolor='white', hatch='x')
    timeout = mpatches.Patch(facecolor='#CBCBCB', label='Timeout', edgecolor='white', hatch='//')
    plt.legend(handles=[complete, failed, timeout], loc='right', bbox_to_anchor=(1.27, 0.8))
    ax[0].set_ylabel('Run Success\n', fontsize=26)
    plt.tight_layout()



def plot_evaluation_metrics(data):
    groups = {'low': 'Low', 'med': 'Medium', 'high': 'High'}
    output = {'Clust_ARI': 'Mutation\nClustering ARI', 'PresAbs_Precision': 'Mutation\nPresence Precision', 'Tree_Accuracy': 'Ancestral Relationship\nAccuracy'}
    plt.rcParams.update({'font.size': 20, 'xtick.bottom': False, 'ytick.left': False})
    fig, ax = plt.subplots(len(output),len(groups), sharex=False, figsize=(12*len(groups),5*(len(output))))
    for i in range(len(output)):
        for j in range(len(groups)):
            if list(output.keys())[i] != 'fail':
                sns.boxplot(y='Method', x=list(output.keys())[i], data=data[data.Sample_groups == list(groups.keys())[j]], ax=ax[i,j], order=['TRACERx', 'Pyclone', 'LICHeE', 'CITUP', 'PhyloWGS'], linewidth=2)
                sns.stripplot(y='Method', x=list(output.keys())[i], data=data[data.Sample_groups == list(groups.keys())[j]], ax=ax[i,j], order=['TRACERx', 'Pyclone', 'LICHeE', 'CITUP', 'PhyloWGS'], edgecolor='black', linewidth=1, size=12, alpha=0.5, palette='tab10')
                ax[i,j].set_xlim(0,1.01)
                if i != len(output)-1:
                    ax[i,j].set(xticklabels=[], xlabel='')
                else:
                    ax[i,j].set_xlabel('Performance', fontsize=26)
            if i == 0:
                ax[i,j].set_title('Sample Group: {}'.format(groups[list(groups.keys())[j]]), fontdict={'fontsize':28})
            if j != 0:
                ax[i,j].set_ylabel('')
            else:
                ax[i,j].set_ylabel(output[list(output.keys())[i]], fontsize=26)

            if j == 0:
                if list(output.keys())[i] in ['Clust_ARI', 'mut_clust_ARI'] or 'PresAbs' in list(output.keys())[i]:
                    ax[i,j].set_yticks(ax[i,j].get_yticks(), labels=['TRACERx', 'Pyclone', 'LICHeE', 'CITUP', 'PhyloWGS'])
                else:
                    ax[i,j].set_yticks(ax[i,j].get_yticks(), labels=['TRACERx', 'Pyclone +\nTRACERx', 'LICHeE', 'CITUP', 'PhyloWGS'])
            else:
                ax[i,j].set_yticklabels([])
                
    plt.tight_layout()
