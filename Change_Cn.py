import MDAnalysis as mda
from MDAnalysis.lib import NeighborSearch, distances
import pandas as pd
import numpy as np
import warnings
import time
import os
import streamlit as st
warnings.filterwarnings("ignore")


@st.cache_resource()
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv(index=False).encode('utf-8')


@st.cache_resource()
def read_file(file):
    u = mda.Universe(file, type='xyz')
    return u

def change_numer(data, frame, save_name):
    data1 = pd.read_csv(data, index_col='Atom Index')
    num = len(data1)
    df = pd.DataFrame(columns=['帧数', 'unchanged_number', 'change_numer'])
    for i in range(0, frame):
        count = 0
        for j in range(num):
            first_num = data1.iloc[j, 0]
            now_num = data1.iloc[j, i]
            if now_num == first_num:
                count += 1
        change_num = num - count
        df.loc[i] = [i + 1, count, change_num]
        df.to_csv(save_name, index=False, encoding='utf-8-sig')


def save_to_csv(df, neighbor_df, atom_indices, save_name, neighbor_save_name):
    # 将原子序号添加到 DataFrame 的第一列
    df.insert(0, 'Atom Index', atom_indices + 1)
    # 将 DataFrame 对象保存为 CSV 文件
    df.to_csv(save_name, index=False)
    # 将邻居氧原子序号保存到相应的数据框列中
    neighbor_df.insert(0, 'Atom Index', atom_indices + 1)
    neighbor_df.to_csv(neighbor_save_name, index=False)


st.title('配位数变化')
# 选择原子
cen_atom = st.selectbox(
    '中心原子选择',
    ('Si', 'Al', 'B', 'Na', 'Mg', 'Ca', 'K', 'Li', 'Se', 'Ba', 'O'))
st.write('你选择:', cen_atom)
uploaded_file = st.file_uploader("Choose a file")
if uploaded_file is not None:
    dataframe = pd.read_csv(uploaded_file)
    st.write(dataframe)

if st.button('统计变化信息'):
    progress_text = "Operation in progress. Please wait."
    my_bar = st.progress(0, text=progress_text)
    previous_frame = dataframe.iloc[:, 1]
    df = pd.DataFrame(columns=['frame', 'decrease', 'stable', 'increase'])
    for i in range(1, len(dataframe.columns)):
        now_frame = dataframe.iloc[:, i]
        decrease_count = 0
        unchanged_count = 0
        increase_count = 0
        if i == 0:
            pass
        else:
            for j in range(len(now_frame)):
                diff = now_frame[j] - previous_frame[j]
                if diff < 0:
                    decrease_count += 1
                elif diff > 0:
                    increase_count += 1
                else:
                    unchanged_count += 1
        df.loc[i] = [i, decrease_count, unchanged_count, increase_count]
        my_bar.progress((i-1) / (len(dataframe.columns)-2), text=progress_text)
    time.sleep(1)
    coordination_csv = convert_df(df)
    st.download_button(
        label="Download coordination as CSV",
        data=coordination_csv,
        file_name=f'统计结果_{cen_atom}.csv',
        mime='text/csv',
    )
else:
    st.write('Goodbye')