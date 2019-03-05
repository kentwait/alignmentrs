import re


sp_sub = re.compile(r'\s+')
col_meta_regexp = re.compile(r'meta\|(\S+)\=(\S+)')


def make_col_meta_string(column_metadata, included_keys, encoders):
    included_values = (
        (k, v) for k, v in column_metadata.to_dict(orient='list').items()
        if k in included_keys
    )
    return ' '.join([
        col_meta_tostr(k, v, encoders[k] if k in encoders.keys() else str)
        for k, v in included_values
    ])

def col_meta_tostr(key, value, encoder: callable):
    return 'meta|{}={}'.format(key, sp_sub.sub('', encoder(value)))

def make_col_meta_dict(description, decoders):
    matches = col_meta_regexp.findall(description)
    return {
        k: (decoders[k](v) if k in decoders.keys() else eval(v))
        for k, v in matches
    }
