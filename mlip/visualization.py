import pandas as pd
import altair as alt
    

def replicate_variant_plot(input_tsv_filepath, output_html_filepath):
        df = pd.read_csv(input_tsv_filepath, sep='\t')
        chart = alt.Chart(df, width=400, height=400).mark_point().encode(
            x=alt.X('Frequency_1:Q', title='Sample 1 variant frequency'),
            y=alt.Y('Frequency_2:Q', title='Sample 2 variant frequency'),
            tooltip=[
                alt.Tooltip('segment:N', title='Segment'),
                alt.Tooltip('position:Q', title='Position'),
                alt.Tooltip('alt:N', title='Alternative'),
                alt.Tooltip('Frequency_1:Q', title='Frequency 1'),
                alt.Tooltip('Frequency_2:Q', title='Frequency 2')
            ]
        ).interactive().properties(
            title=alt.TitleParams(
                text='Replicate variant calls',
                fontSize=20,
                anchor='middle'
            ),
            padding={'left': 20, 'right': 20, 'top': 10, 'bottom': 10}
        ).configure_axis(
            labelFontSize=18,
            titleFontSize=18
        )
        chart.save(output_html_filepath)

