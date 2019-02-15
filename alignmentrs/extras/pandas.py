import pandas as pd
from libalignmentrs.position import Block, BlockSpace
from alignmentrs.aln.classes import Alignment, CatAlignment


def blockspace_to_df(aln: Alignment):
    block_space = aln._linspace
    df_dict = {'start': [], 'stop': [], 'state': []}
    for block in block_space.to_blocks():
        df_dict['start'].append(block.start)
        df_dict['stop'].append(block.stop)
        df_dict['state'].append(block.id)

    return pd.DataFrame(df_dict)


def cat_blockspace_to_df(cat_aln: CatAlignment):
    cat_space = cat_aln._linspace
    sub_spaces = cat_aln._subspaces

    df_dict = {
        'cat_start': [], 'cat_stop': [],
        'aln_start': [], 'aln_stop': [],
        'aln_name': [], 'state': [],
    }
    for catblock in cat_space.to_blocks():
        for block in sub_spaces[catblock.id].to_blocks():
            df_dict['cat_start'].append(catblock.start)
            df_dict['cat_stop'].append(catblock.stop)
            df_dict['aln_name'].append(catblock.id)
            df_dict['aln_start'].append(block.start)
            df_dict['aln_stop'].append(block.stop)
            df_dict['state'].append(block.id)
    
    return pd.DataFrame(df_dict)

