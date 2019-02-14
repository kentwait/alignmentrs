import pandas as pd
from alignmentrs import Block, BlockSpace


def blockspace_to_df(block_space: BlockSpace, state_map=None):
    df_dict = {'start': [], 'stop': [], 'state': []}
    for block in block_space.to_blocks():
        df_dict['start'].append(block.start)
        df_dict['stop'].append(block.stop)

        state = block.id if state_map is None else state_map[block.id]
        df_dict['state'].append(state)

    return pd.DataFrame(df_dict)
