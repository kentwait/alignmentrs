import re


sp_sub = re.compile(r'\s+')

def make_col_meta_string(column_metadata, included_keys, encoders: dict=None):
    if encoders is None:
        encoders = {}
    included_values = (
        (k, v) for k, v in column_metadata.to_dict(orient='list').items()
        if k in included_keys
    )
    return ' '.join([
        col_meta_strfmt(k, v, encoders[k] if k in encoders.keys() else str)
        for k, v in included_values
    ])

def col_meta_strfmt(key, value, encoder: dict):
    return 'meta|{}={}'.format(key, sp_sub.sub('', encoder(value)))
