import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def plot_mcc_density (sample_1, save_plot = False, save_name = 'CC_mmc.png'):
    for sample in sample_1:
        df = pd.read_csv(f'D:\Jupyter_Notebook\Xenium_Jupyter_notebook\Correlation_Mapping\{sample}_CorrelationMapping.csv', comment = "#")
        df['sample'] = df['cell_id'].map(lambda name: name.split('_')[0])
        df_temp = df.filter(['sample', 'subclass_correlation_coefficient'])
        if "df_all" not in locals():
            df_all = df_temp
        else:
            df_all = pd.concat([df_all, df_temp])

    dict_ = dict(zip(df_all['sample'],df_all['sample']))
    dict_.update({"Region1": "february-test","S1" : "march-test"})
    df_all['sample'] = df_all['sample'].map(dict_)

    order_sample = df_all['sample'].unique() ## All samples

    ax = sns.histplot(data = df_all, x = 'subclass_correlation_coefficient', hue='sample', hue_order= order_sample,
             element="step", cumulative= True, fill= False, common_norm=False,
             stat='density')
    plt.xlim(0,0.85)
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))

    if save_name == 'CC_mmc.png':
        save_name = f'CC_{sample_1[0]}.png'

    if save_plot:
        plt.savefig('Gallery/save_name')