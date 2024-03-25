import streamlit as st
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib import distances
import math
import pandas as pd
import time
import os


@st.cache_resource()
def convert_df(df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return df.to_csv().encode('utf-8')


@st.cache_resource()
def read_file(file):
    u = mda.Universe(file, type='xyz')
    return u


st.title('键长分析')
# 选择原子
cen_atom = st.selectbox(
    '中心原子选择',
    ('Si', 'Al', 'B', 'Na', 'Mg', 'Ca', 'K', 'Li', 'Se', 'Ba', 'O'))
st.write('你选择:', cen_atom)
sur_atom = st.selectbox(
    '周围原子选择',
    ('Si', 'Al', 'B', 'Na', 'Mg', 'Ca', 'K', 'Li', 'Se', 'Ba', 'O' ))
st.write('你选择:', sur_atom)
cutoff = st.number_input('请输入截断半径', step=0.001)
st.write('截断半径：', cutoff)

# 设置一个临时目录用于保存上传的文件
TEMP_DIR = "temp_files"
os.makedirs(TEMP_DIR, exist_ok=True)

# 获取上传的文件并保存到临时目录中
uploaded_files = st.file_uploader("Choose a xyz file", accept_multiple_files=True, type='xyz')
uploaded_file_paths = []  # 存储上传文件的路径
for uploaded_file in uploaded_files:
    file_path = os.path.join(TEMP_DIR, uploaded_file.name)
    with open(file_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    uploaded_file_paths.append(file_path)  # 将文件路径添加到列表中

# 在临时目录中找到上传的文件，并显示文件路径
for file_path in uploaded_file_paths:
    st.write("Uploaded file path:", file_path)

if st.button('统计键长信息'):
    bond_length_data = []
    progress_text = "Operation in progress. Please wait."
    my_bar = st.progress(0, text=progress_text)
    for index, uploaded_file in enumerate(uploaded_file_paths):
        time.sleep(0.01)
        u = read_file(uploaded_file)
        atom1_atoms = u.select_atoms(f'name {cen_atom}')
        atom2_atoms = u.select_atoms(f'name {sur_atom}')
        bond_length = {}
        for atom1_atom in atom1_atoms:
            atoms_distances = distances.distance_array(atom1_atom.position, atom2_atoms.positions)
            indices = np.where(atoms_distances <= cutoff)[1]
            atom2_atom = atom2_atoms[indices]
            for atom2 in atom2_atom:
                bond_lengthes = atom1_atom.position - atom2.position
                bond_mid = []
                for i in bond_lengthes:
                    x = i ** 2
                    bond_mid.append(x)
                length = math.sqrt(sum(bond_mid))
                bond_length[(atom1_atom.index + 1, atom2.index + 1)] = length
        my_bar.progress((index + 1) / len(uploaded_file_paths), text=progress_text)
        length_series = pd.Series(bond_length)
        length_series.name = f'{uploaded_file}'
        bond_length_data.append(length_series)
    time.sleep(1)
    df_bond_length = pd.concat(bond_length_data, axis=1)
    bond_length_csv = convert_df(df_bond_length)
    st.download_button(
        label="Download data as CSV",
        data=bond_length_csv,
        file_name=f'{cen_atom}_bond_length_data.csv',
        mime='text/csv',
    )
else:
    pass


# st.write('致谢：本脚本由2023级研究生康泽玉，2020级本科生李新宇开发')

