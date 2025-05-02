'''
    Filename: plot_passage_pie.py
    Notes: Makes pie chart with PASSAGE filter/orientation stats, using uinfo from the spreadsheet
    Author : Ayan
    Created: 02-05-25
    Example: run plot_passage_pie.py
'''
from header import *
from util import *

start_time = datetime.now()

# -------------------------------------------------------------------------------------------------------
def plot_sunburst_pie(df, args, outer_col='nPA', inner_col='filters'):
    '''
    To plot concentric pie chart given a df, for a bunch of criteria
    Plots and saves the figure
    '''
    df['count'] = 1

    fig = px.sunburst(df, path=[inner_col, outer_col], values='count', hover_data=['count'])
    fig.update_traces(textinfo='label+value+percent entry', texttemplate='<b>%{label}</b><br>(%{value})<br>%{percentEntry:.1%}')

    if args.fortalk: fig.update_layout({'plot_bgcolor': 'rgba(0, 0, 0, 0)', 'paper_bgcolor': 'rgba(0, 0, 0, 0)'})
    df = df.drop('count', axis=1)
   
    # ----------annotate and save the diagram----------
    figname = args.output_dir / 'plots' / f'passage_allfields_sunburst_diagram.pdf'
    fig.write_image(figname)
    print(f'\nSaved figure as {figname}')
    #fig.show()

    return fig

# -------------------------------------------------------------------------------------------------------
def read_passage_spreadsheet(args=None, filename=None):
    '''
    Reads PASSAGE spreadsheet and and returns a dataframe
    '''
    if filename is None:
        filename =  args.input_dir / 'JWST PASSAGE Cycle 1 - Cy1 Executed.csv'

    df = pd.read_csv(filename)
    df = df[['Par#', 'Obs Date', 'Filter', 'Mode', 'Total Exposure Time [sec]']]
    df = df.fillna(method='ffill')

    df = df[~(df['Par#'] == 'SKIP')]
    df = df[df['Par#'].map(lambda x: int(x[3:])) < 600] # to get only Cy 1 objects

    df = df[~(df['Obs Date'].str.contains('SKIPPED'))]

    return df

# --------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    args = parse_args()
    if not args.keep: plt.close('all')

    # -------reading in spreadsheet-------
    df = read_passage_spreadsheet(args)

    # -----determining filters and orientations-------
    df = df[~(df['Mode'].str.contains('DIRECT'))][['Par#', 'Filter', 'Mode']]
    df = df.drop_duplicates(['Par#', 'Filter', 'Mode'])
    df = df.groupby(['Par#', 'Filter']).agg('count').reset_index()
    df['Mode'] = df['Mode'].apply(lambda x: f'{x} PA' if x <= 1 else f'{x} PAs')
    df = df.pivot_table(index=['Par#', 'Mode'], values=['Filter'], aggfunc=lambda x: ','.join(x)).reset_index()

    # -----plotting sunburst chart-------
    fig = plot_sunburst_pie(df, args, outer_col='Mode', inner_col='Filter')

    print(f'Completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
