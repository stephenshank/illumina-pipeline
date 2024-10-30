import http.server
import socketserver
import webbrowser
import socket

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



def find_free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('', 0))
        return s.getsockname()[1]

def start_server():
    PORT = find_free_port()
    handler = http.server.SimpleHTTPRequestHandler

    with socketserver.TCPServer(("", PORT), handler) as httpd:
        print(f"Serving on port {PORT}")
        webbrowser.open(f"http://localhost:{PORT}")
        httpd.serve_forever()

if __name__ == "__main__":
    start_server()

