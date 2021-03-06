"""
This file is part of Oligo-ASST, released under the 3-Clause BSD License.
See License.txt or go to https://opensource.org/licenses/BSD-3-Clause for
full details.

Copyright 2020 Miler T. Lee.


application.py: Oligo-ASST web interface using Dash framework
Author: Miler Lee
"""

#import os
import base64
import io

import dash
from dash.dependencies import Input, Output, State
import dash_html_components as html
import dash_core_components as dcc
import dash_table
import dash_bio
#from dash_bio_utils import protein_reader as pr  #dependencies bug
import pandas as pd

import oligo_asst
import mtl_fasta

app = dash.Dash(__name__, assets_folder='assets/')
server = app.server
#application = app.server #for AWS


def about_blurb():
	"""
	Description of the web app	
	"""
	blurb = [html.P('Oligo-ASST (Anti-Sense Spaced Tiled) designs antisense oligos for RNaseH-mediated RNA digestion, for applications such as rRNA depletion of RNA-seq libraries.'),
	html.P(html.Img(src='assets/phelps20_graphical_abstract_v2.png', alt = 'Oligo-ASST diagram', style= {'width': '75%'}), style={'text-align': 'center'}),
	html.P('Given an input sequence provided as a FASTA file, Oligo-ASST attempts to find oligo target regions sparsely tiled across the length of the sequence as close to an RNA-DNA melting temperature of 70-80C as possible.'),
	html.P(html.Img(src='assets/oligo_parameter_cartoon.png', alt = 'Oligo parameter diagram', style= {'width': '75%'}), style={'text-align': 'center'}),	
	html.P(['If multiple input sequences are provided in the FASTA file, Oligo-ASST can optionally design oligos that target multiple sequences, to decrease the total number of oligos needed for full targeting. This feature works best if the sequences are not highly divergent and have been aligned by a program such as ', html.A('MUSCLE', href='https://www.ebi.ac.uk/Tools/msa/muscle/'), ' beforehand.']),
	html.P('If you find this useful, please cite Phelps et al, 2020. Optimized design of antisense oligomers for targeted rRNA depletion, bioRxiv.'),
	html.P(['Contact: ', html.A('M.T. Lee Lab', href='https://mtleelab.pitt.edu'), ', miler at pitt.edu.'])
	]


	return blurb


def precalc_oligos():
	"""
	Populates the precal oligo panel
	"""
	
	return html.Div(children = [
		html.Table([
			html.Tr([
				html.Td(["Xenopus laevis"]),
				html.Td(
					[html.A(href= '/assets/xla_oligos_details.csv',
					children='[oligos+details]',
					download='xlaevis_rrna_oligos_details.csv'),
					
					" (Verified, Phelps 2020)"					
				] ,style = {'padding': '10px'}),
					
			], style = {'border-bottom': '1px solid gray'}),

			html.Tr([
				html.Td(["Xenopus tropicalis"]),
				html.Td([
					html.A(href= '/assets/xtr_oligos.csv',
					children='[oligos]',
					download='xtropicalis_rrna_oligos.csv'),
					
					html.A(href= '/assets/xtr_full_details.csv',
					children='[full details]',
					download='xtropicalis_rrna_full_details.csv'),
					
					" (Predicted)"
				] ,style = {'padding': '10px'}),
			], style = {'border-bottom': '1px solid gray'}),
			
			
			html.Tr([
				html.Td(["Zebrafish (maternal+somatic)"]),
				html.Td(
					[html.A(href= '/assets/zf_oligos_details.csv',
					children='[oligos+details]',
					download='drerio_rrna_oligos_details.csv'),
					
					" (Verified, Phelps 2020)"
					
				],style = {'padding': '10px'})				
			], style = {'border-bottom': '1px solid gray'}),
		], style = {'border-collapse': 'collapse'})
			
	])


def table_component():
	"""
	Returns the children for a clean Html Div for the oligo tables
	"""

	return [
			html.Div(id='table-buttons',
				style= {'display': 'none'},
				children=[
						html.Table(children = html.Tr(
							children=[
							html.Td(children=html.H4(id='oligo-count', children=[])),
							html.Td(dcc.RadioItems(
								id='table-display-radio',
								options = [
									{'label': 'Oligo list', 'value': 'oligo-summary'},
									{'label': 'Full details', 'value': 'oligo-details'}
								],
								value='oligo-summary',
								labelStyle={'display': 'inline-block'}
								),
							),
							html.Td(children=html.Abbr('\u24D8', title='"Oligo list" contains all the unique antisense oligos across all targeted sequences. "Full details" lists all targeted positions, so oligos may appear more than once.', style={'display': 'inline-block'}))
						]))
				]

			),
		
			html.Div(id='table-container',
					children = [
					],
					style={'width': '100%', 'display': 'inline-block', 'vertical-align': 'top', 'padding-left': 10}
			)
		]


def format_oligo_detail_table(a, display = True):
	"""
	Returns a dash_table
	"""

	dt = dash_table.DataTable(id='oligo-detail-table',
		columns= [{"name": i, "id": i} for i in a.columns],
		
		data=a.to_dict('records'),
		style_header={'backgroundColor': 'rgb(230,230,230)', 'fontWeight':'bold'},

		style_cell={'width': '70px'},
		
		style_cell_conditional=[
			{'if': {'column_id': 'Name'}, 'textAlign': 'left'},
			{'if': {'column_id': 'Target'}, 'textAlign': 'left'},			
			{'if': {'column_id': 'Antisense_oligo'}, 'textAlign': 'left'},
			{'if': {'column_id': 'Target_seq'}, 'textAlign': 'left'}
		],

		style_data_conditional=[
			{'if': {'row_index': 'odd'}, 'backgroundColor': 'rgb(248,248,248)'},
			{'if': {'column_id': 'Tm', 'filter_query': '{Tm} lt 65'}, 'color': 'brown', 'fontWeight':'bold', 'backgroundColor': 'yellow'},		
			{'if': {'column_id': 'Tm', 'filter_query': '{Tm} lt 55'}, 'color': 'white', 'fontWeight':'bold', 'backgroundColor': 'orange'},
		],
		
#		style_table={'maxHeight': '800px', 'overflowY': 'scroll'},
		style_table={'overflowY': 'scroll'},
		fixed_rows={ 'headers': True, 'data': 0 },
		
		sort_action='native',
		export_format='csv',
		export_headers='names',
	)
	
	if display:
		display_string = 'block'
	else:
		display_string = 'none'

	return html.Div(id='table-full-container', style={'display': display_string}, children=dt)


def format_oligo_summary_table(a, display = True):
	"""
	Returns a dash_table with only the columns for Name, oligo, and wildcard
	"""
	
	summary = a.filter(items = ['Name', 'Antisense_oligo', 'Tm', 'N_targets', 'Wcards'])
	summary = summary.sort_values(by=['Wcards', 'Name']).drop_duplicates()
	num_oligos = len(summary)

	dt = dash_table.DataTable(id='oligo-summary-table',
		columns= [{"name": i, "id": i} for i in summary.columns],
		
		data=summary.to_dict('records'),
		style_header={'backgroundColor': 'rgb(230,230,230)', 'fontWeight':'bold'},

		style_cell={'width': '60px', 'textAlign': 'left'},
		
		style_data_conditional=[
			{'if': {'row_index': 'odd'}, 'backgroundColor': 'rgb(248,248,248)'},
			{'if': {'column_id': 'Tm', 'filter_query': '{Tm} lt 65'}, 'color': 'brown', 'fontWeight':'bold', 'backgroundColor': 'yellow'},			
			{'if': {'column_id': 'Tm', 'filter_query': '{Tm} lt 55'}, 'color': 'white', 'fontWeight':'bold', 'backgroundColor': 'orange'},
		],		
		
		style_table={'overflowY': 'scroll'},
		fixed_rows={ 'headers': True, 'data': 0 },
		
		sort_action='native',
		export_format='csv',
		export_headers='names',
	)
	
	if display:
		display_string = 'block'
	else:
		display_string = 'none'

	return html.Div(id='table-summary-container', style={'display': display_string}, children=dt), num_oligos



def sequence_index_map(gapped_sequence):
	"""
	Creates a list that maps index positions in the non-gapped
	version of a sequence to the positions in the gapped sequence
	"""
	return [i for i,x in enumerate(gapped_sequence+'$') if x != '-']


def set_coverage(ident, seq, data):

	"""
	Based on DashTable, sets coverage for the sequence
	Assumes sorted. Prevents overlaps
	"""
	cov = []
	last_end = 0
	
	index_map = sequence_index_map(seq)

	for row in sorted(data, key = lambda x: int(x['Start'])):
		if row['Target'] != ident:
			continue
	
		tooltip = '%d..%d %s ' % (row['Start'], row['End'], row['Name'])
		if 65 <= row['Tm']:
			bgcolor = 'lightgreen'

		elif 55 <= row['Tm'] < 65:
			tooltip += '(Tm between 55-65)'
			bgcolor = 'yellow'

		else:
			tooltip += '(Tm < 55)'
			bgcolor = 'orange'

		if(int(row['Start']) <= int(row['End'])):
			start_coord = index_map[max(last_end, int(row['Start'])-1)]
			end_coord = index_map[max(last_end, int(row['End']))]
		
			cov.append({'start': start_coord, 'end': end_coord, 'bgcolor': bgcolor, 'tooltip': tooltip, 'underscore': True})
			last_end = max(int(row['End']), last_end)

	return cov


###
# Layout components
###

def upload_alert_blurb(fname = None, msg = None):
	"""
	"""
	if not fname:
		return ['No file loaded. ', html.A('Load a sample file', href='#', id='load-sample-fasta-link')]
	elif msg:
		return ['{}: {} '.format(msg, fname), html.A('Load a sample file', href='#', id='load-sample-fasta-link')]
	else:
		return ['{} uploaded (Click Sequence tab to view)'.format(fname), html.A('', href='#', id='load-sample-fasta-link')]


def layout_fasta_upload_div():
	return html.Div(
		id='seq-view-fasta-upload',
		children=[
			dcc.Upload(
				id='upload-fasta-data',
				max_size=1049000,
				style={
						'width': '95%',
						'height': '75px',
						'lineHeight': '75px',
						'borderWidth': '1px',
						'borderStyle': 'solid',
						'borderRadius': '5px',
						'textAlign': 'center',
						'margin': '5px'
				},
				children=html.Div([
					"Drag and drop or click to upload a \
					FASTA file (<1MB)"
				]),
			),
			html.Div(id='upload-alert', children=upload_alert_blurb(), style={'font-size': '0.8em', 'font-style':'italic','display':'inline-block'})

		]
	)


def layout_basic_params_div():
	return html.Div([
	
		html.P('Default: 39-40nt oligos spaced <=30 nts apart, favoring melting temperatures between 70-80C. Some parameter combinations may not be feasible given the input sequence', style = {'margin-top': '25px'}),
		
		html.Table(children=[
			html.Tr([
				html.Th('Min oligo length'),
				html.Th('Max oligo length'),
				html.Th('Max untiled length')									
			]),
			html.Tr([
				html.Td(dcc.Input(
					id='min-oligo-input',
					type='number',
					min=5,
					max=1000,
					value='39'),
				),
				html.Td(dcc.Input(
					id='max-oligo-input',									
					type='number',
					min=5,
					max=1000,
					value='40')
				),
				html.Td(dcc.Input(
					id='max-untiled-input',
					type='number',
					min=5,
					max=1000,
					value='30')
				)
			]),
			html.Tr([
				html.Td(html.H4('Tm range (lo/hi C)'),
				),
				html.Td([dcc.Input(
					id='tm-low-input',
					type='number',
					min=40,
					max=95,
					value='70'
					),
				]),
				html.Td([dcc.Input(
					id='tm-high-input',
					type='number',
					min=40,
					max=95,
					value='80'
					),
				]),
			]),
			
			html.Tr([
				html.Td(html.H4('Oligo name prefix'),
				),
				html.Td([dcc.Input(
					id='oligo-prefix-input',
					type='text',
					placeholder='oligo',
					value='oligo',
					debounce=True
					),
#					html.Abbr('?', title='Oligos will be named consecutively as "<PREFIX>_<NUM>"')
				]),
				html.Td(
#					html.Button(id='go-button2', children=html.B('Calculate')),
#					style = {'text-align': 'right'}
				)
			])
		]),	
	])
	

def layout_multiseq_params_div():
	return html.Div(children=[

		html.Table([
			html.Tr([			
				html.Th(['For >1 input sequence ', html.Abbr('\u24D8', title='For best performance when calculating shared oligos, align your sequences using a tool such as CLUSTAL or MUSCLE and output as a FASTA file prior to uploading here')],
						colSpan=2,
						style = {'textAlign': 'left'}),
				html.Th()
			]),
			
			html.Tr([
				html.Td(dcc.RadioItems(
						id='multi-mode-radio',
						options = [
							{'label': 'Independently design oligos for each sequence', 'value': 'basic'},
							{'label': 'Sequences are not aligned, design shared oligos', 'value': 'combined'},
							{'label': 'Sequences are aligned, design shared oligos', 'value': 'aligned'},
							
						],
						value='basic',
						labelStyle={'display': 'inline-block'}
						)
				, colSpan=3)
			]),
		]),

		html.Div(id='shared-oligo-params', children=[
			html.Table([			
				html.Tr([
					html.Td([], style={'width': '30px'}),
					html.Td([html.Div('Max wildcard expansions', style={'font-weight': 'normal'})]),
					html.Td([
						dcc.Input(
							id='num-wildcards-input',
							type='number',
							min=1,
							max=16,
							placeholder='1 = no wildcards'
						),
						html.Abbr("\u24D8", title='Wildcards represent combinations of bases (e.g., R means A or G, N means any base). Expansions are the number of targets per oligo, as a function of number of wildcards encoded. E.g., 0 wildcards = 1 expansion, 2 wildcards each representing 2 possible bases is 2x2 = 4 total expansions')
					])
				]),
				
				html.Tr([
					html.Td([]),
					html.Td([
						dcc.Checklist(
							id='trim-wildcards-checkbox',
							options=[{'label': 'Trim terminal wildcards','value': 'trim'}],
							value=['trim']
						)
					]),
				
					html.Td([html.Abbr('\u24D8', title='If an oligo begins or ends with a wildcard base, trim it off even if it violates the oligo length range')
					])
				])
			])
		])
		
	])
					


def layout():  
	return html.Div([

	#Left side: Tabbed: FASTA entry / sequence viewer
	html.Div(id='tabs',
#		className='column',
		style= {'width': '40%', 'display': 'inline-block', 'vertical-align' : 'top'},
		children=[
		html.H2('Oligo-ASST', style={'display':'inline-block'}),
		html.H5('Design antisense oligos for RNA depletion',style={'display':'inline', 'margin-left': '10px', 'color': '#758b9e'}),
		dcc.Tabs(id='seq-view-tabs', value='ctrl-tab', 
			children=[
			#FASTA entry tab
			dcc.Tab(
				label='Parameters',
				value='ctrl-tab',
				id='ctrl-tab',
				style={'padding': '10px'},
				selected_style={'padding': '10px'},				
				children=html.Div(children=[
					#Fasta upload section
					layout_fasta_upload_div(),
					layout_basic_params_div(),
					layout_multiseq_params_div(),
					
					html.Div(
						[html.P(''),
						html.Button(id='go-button', children=html.B('Calculate')),
						html.I(id='go-status-msg', children='')					  
						],
					),
				])
			),
		
			#Sequence viewer tab
			dcc.Tab(
				label='Sequence',
				value='seq-tab',
				style={'padding': '10px'},
				selected_style={'padding': '10px'},								
				children=html.Div(id='seq-view-container',
					children=[

						##Selector (top)
						html.Div(id='seq-view-info-container',
							children=html.Div(id='seq-view-info',
								children=[html.Div(id='seq-view-info-desc',
									children=[html.Span("",
										className='seq-view-info-element-title'),
										html.Div(id='desc-info', children=[])
									]),
								],

							)
						),
					
						##Sequence viewer (bottom)
						html.Div(id='seq-view-component-container-outer',
							children = html.Div(id='seq-view-component-container', children=[])
						),
						
					],
				),
			),

			#Precalc oligos tab
			dcc.Tab(
				label='Precalculated oligos',
				value='precalc-tab',
				style={'padding': '10px'},
				selected_style={'padding': '10px'},								
				children=precalc_oligos()
			),
			
			#About tab
			dcc.Tab(
				label='About',
				value='about-tab',
				style={'padding': '10px'},
				selected_style={'padding': '10px'},				
				
				children=html.Div(about_blurb(), style={'font-size': '0.9em'})
			)
		]),
		html.Hr(style = {'border-top': '1px','margin-top': '50px'}),
		html.P(['Oligo-ASST is brought to you by the ', html.A('M.T. Lee lab', href='https://mtleelab.pitt.edu'), ', University of Pittsburgh. ', 'If you use this tool, please cite Phelps et al, 2020, Optimized design of antisense oligomers for targeted rRNA depletion, bioRxiv.'], style={'font-size': '0.75em'})
	]),


	#Spacer
	html.Div(id='spacer',
#	className='column',
	style={'width': '2%', 'display': 'inline-block'},
	children=[]
	),


	#Right side: table of oligos
	html.Div(id='right-side',
		style= {'width': '55%', 'display': 'inline-block'},
#		className='column',
		children=table_component()
	),
	dcc.Loading(id='loading-widget',
                  children=[html.Div(id="loading-output")],
                type="default",
                fullscreen = True
    )
	
	])



#		children = ['{} uploaded (Click Sequence tab to view)'.format(fname), html.A('', href='#', id='load-sample-fasta-link')]
			
#	return html.Div(id='upload-alert', children, style={'font-size': '0.8em', 'font-style':'italic','display':'inline-block'})


#			html.Div(id='upload-alert', upload_alert_blurb(), style={'font-size': '0.8em', 'font-style':'italic','display':'inline-block'})



def callbacks(_app):

	#FASTA file upload
	@_app.callback(
		[Output('upload-alert', 'children'),
		Output('seq-view-component-container-outer', 'children'),
		Output('right-side', 'children'),
		],
		[Input('upload-fasta-data', 'filename'),
		Input('upload-fasta-data', 'contents'),
		Input('load-sample-fasta-link', 'n_clicks')]
	)
	def update_upload_msg(fname, contents, sample_file_clicks):
		try:
			if not fname:
				if sample_file_clicks:
					fname = 'zebrafish_28s_aligned.fa'
					content_str = open('assets/' + fname).read()
				else:
					return upload_alert_blurb(), None, table_component()	
			else:
				content_type, content_str = contents.split(',')
				content_str = base64.b64decode(content_str).decode('UTF-8')
		
			seq_viewers = []
			multi_fasta = mtl_fasta.parse_fasta(content_str)
#			multi_fasta = pr.read_fasta(datapath_or_datastring=content_str, is_datafile=False)
			ident_hash = {}
			for i, fasta in enumerate(multi_fasta):
#				sequence = fasta.get('sequence', '')
#				ident = fasta.get('description', {}).get('Header', 'Unnamed').split()[0]
				ident = fasta[0].split()[0] #First word only
				sequence = fasta[1]
				if ident not in ident_hash:
					ident_hash[ident] = 1
				else:
					ident_hash[ident] += 1
					ident = ident + '(' + str(ident_hash[ident]) + ')'
				
				sv = dash_bio.SequenceViewer(id='sequence-viewer-%d' % i,
								sequence = sequence,
								title = ident,
								charsPerLine = 50, search=False,
								badge=False, selection=False,
								coverage=[])				
				seq_viewers.append(sv)
			
			#TODO: Perform a check for proper DNA sequence
			#TODO: Filesize warning?
			container = html.Div(id='seq-view-component-container', children=seq_viewers)

			return upload_alert_blurb(fname), container, table_component(),	
		except Exception as e:
			print(e)
			return upload_alert_blurb(fname, 'File upload error'), None, table_component()
			


	#Do the computation upon button press
	@_app.callback(
		[Output('table-container', 'children'),
		Output('table-buttons', 'style'),
		Output('oligo-count', 'children'),
		Output('loading-output', 'children')
		],
		[Input('go-button', 'n_clicks'),
#		Input('go-button2', 'n_clicks'),
		],
		[State('seq-view-component-container', 'children'),
		State('min-oligo-input', 'value'),
		State('max-oligo-input', 'value'),
		State('max-untiled-input', 'value'),
		State('tm-low-input', 'value'),
		State('tm-high-input', 'value'),		
		State('multi-mode-radio', 'value'),
		State('num-wildcards-input', 'value'),
		State('table-display-radio', 'value'),
		State('oligo-prefix-input', 'value'),
		State('trim-wildcards-checkbox', 'value'),
#		State('aligned-checkbox', 'value')
		]
	)

	def do_tiling(n_clicks, seq_viewer, min_len_str, max_len_str, max_untiled_len, tm_low, tm_high, multi_mode, n_wildcards, display_mode, oligo_name_prefix, trim_terminal_wildcards):
	
		if not seq_viewer:
			return None, {'display': 'none'}, None, True
		if not n_wildcards:
			n_wildcards = 1
		if n_clicks:
			min_len = min(int(min_len_str), int(max_len_str))
			max_len = max(int(min_len_str), int(max_len_str))
			min_tm = min(int(tm_low), int(tm_high))
			max_tm = max(int(tm_low), int(tm_high))
			
			idents = [x['props']['title'] for x in seq_viewer]
			seqs = [x['props']['sequence'] for x in seq_viewer]
			oligos = oligo_asst.do_all(idents, seqs, min_len = min_len, max_len = max_len, tm_low = min_tm, tm_high = max_tm, max_untiled_len = int(max_untiled_len), do_consensus = (multi_mode == 'combined' or multi_mode == 'aligned'), max_wildcard_expansions = int(n_wildcards), oligo_name_prefix = oligo_name_prefix, trim_terminal_wildcards = len(trim_terminal_wildcards), seqs_aligned = (multi_mode == 'aligned'))
			
			try:
				df = pd.read_csv(io.StringIO(oligos), sep=",", dtype = {'Start': 'Int64', 'End': 'Int64'})				
				detail_table = format_oligo_detail_table(df, display_mode == 'oligo-details')
				oligo_table, num_oligos = format_oligo_summary_table(df, display_mode == 'oligo-summary')
				
				return [detail_table, oligo_table], {'display': 'block'}, '%d oligos' % num_oligos, True
			except Exception as e:
				print(e)
				return html.H3('No results -- try different length parameters'), {'display': 'none'}, None, True
			
	#Update sequence highlighting upon calculation
	@_app.callback(
		#Output('sequence-viewer', 'coverage'),
		Output('seq-view-component-container', 'children'),	
		[Input('oligo-detail-table', 'data')],
		[State('seq-view-component-container', 'children')]
	)
	def update_highlights(data_table, seq_viewers):
		new_seq_viewers = []
		for sv in seq_viewers:
			new_sv = dash_bio.SequenceViewer(id=sv['props']['id'],
								sequence = sv['props']['sequence'],
								title = sv['props']['title'],
								charsPerLine = 50, search=False,
								badge=False, selection=False,
								coverage=set_coverage(sv['props']['title'], sv['props']['sequence'], data_table)
						)			
			new_seq_viewers.append(new_sv)
			
		return new_seq_viewers


	@_app.callback(
		[Output('table-summary-container', 'style'),
		Output('table-full-container', 'style')],
		[Input('table-display-radio', 'value')]
	)
	def toggle_table_view(display_mode):
		if display_mode == 'oligo-details':
			return {'display': 'none'}, {'display': 'block'}
		else:
			return {'display': 'block'}, {'display': 'none'}

	@_app.callback(
		Output('shared-oligo-params', 'style'),
		[Input('multi-mode-radio', 'value')]
	)
	def toggle_shared_params_view(multi_mode):
		if multi_mode == 'aligned':
			return {'display': 'block'}
		else:
			return {'display': 'none'}



##### RUNNING
    
app.scripts.config.serve_locally = True
app.config['suppress_callback_exceptions'] = True
app.layout = layout()
app.title = 'Oligo-ASST'
callbacks(app)

					
if __name__ == '__main__':
#	application.run(debug=True, port=8080) #for AWS
	server.run(debug=True)

	
