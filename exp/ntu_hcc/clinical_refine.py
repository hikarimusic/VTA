import pandas as pd
import numpy as np

df = pd.read_csv('clinical.tsv', sep='\t')
df = df.drop_duplicates(subset=['case_id'], keep='first')

gender_mapping = {
    'male': 'male',
    'female': 'female'
}
df.loc[:, 'Gender'] = df['gender'].map(gender_mapping)


race_mapping = {
    'asian': 'asian',
    'white': 'white',
    'black or african american': 'black'
}
df.loc[:, 'Race'] = df['race'].map(race_mapping)


tumor_mapping = {
    'T1': 'T1',
    'T2': 'T2',
    'T2a': 'T2',
    'T2b': 'T2',
    'T3': 'T3',
    'T3a': 'T3',
    'T3b': 'T3',
    'T4': 'T4'
}
df.loc[:, 'Tumor'] = df['ajcc_pathologic_t'].map(tumor_mapping)


nodes_mapping = {
    'N0': 'N0',
    'N1': 'N1',
    'NX': 'NX'
}
df.loc[:, 'Nodes'] = df['ajcc_pathologic_n'].map(nodes_mapping)


metastasis_mapping = {
    'M0': 'M0',
    'M1': 'M1',
    'MX': 'MX'
}
df.loc[:, 'Metastasis'] = df['ajcc_pathologic_m'].map(metastasis_mapping)


stage_mapping = {
    'Stage I': 'Stage I',
    'Stage II': 'Stage II',
    'Stage III': 'Stage III',
    'Stage IIIA': 'Stage III',
    'Stage IIIB': 'Stage III',
    'Stage IIIC': 'Stage III',
    'Stage IV': 'Stage IV',
    'Stage IVA': 'Stage IV',
    'Stage IVB': 'Stage IV'
}
df.loc[:, 'Stage'] = df['ajcc_pathologic_stage'].map(stage_mapping)


event_mapping = {
    'Dead': 1,
    'Alive': 0
}
df.loc[:, 'Event'] = df['vital_status'].map(event_mapping)


df.loc[:, 'Days'] = np.where(
    df['Event'] == 1,
    df['days_to_death'],
    df['days_to_last_follow_up']
)
df.loc[:, 'Days'] = pd.to_numeric(df['Days'], errors='coerce')


cols = df.columns.tolist()
cols[0], cols[1] = cols[1], cols[0]
df = df[cols]
df.iloc[:, 0] = df.iloc[:, 0].str.replace('-', '_')


df.to_csv('LIHC.tsv', sep='\t', index=False)
