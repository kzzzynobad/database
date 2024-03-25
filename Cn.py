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


def cn_o_num(atoms, cutoff, sur_atoms):  # 中心原子 截断半径 被寻原子
    coordination_numbers = []
    neighbor_indices = []
    for atom in atoms:
        ns = NeighborSearch.AtomNeighborSearch(atoms)
        atoms_distances = distances.distance_array(atom.position, sur_atoms.positions)
        distance_cutoff = cutoff
        coordination_number = np.sum(atoms_distances <= distance_cutoff)
        coordination_numbers.append(coordination_number)
        indices = np.where(atoms_distances <= distance_cutoff)[1]
        neighbor_indices.append(list(sur_atoms.indices[indices]))
    return coordination_numbers, neighbor_indices


st.title('配位数')
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


if st.button('统计配位信息'):
    bond_angle_data = []
    progress_text = "Operation in progress. Please wait."
    my_bar = st.progress(0, text=progress_text)
    data = {}
    neighbor = {}
    df = pd.DataFrame(data)
    neighbor_df = pd.DataFrame(neighbor)
    for index, uploaded_file in enumerate(uploaded_file_paths):
        time.sleep(0.01)
        u = read_file(uploaded_file)
        atom1_atoms = u.select_atoms(f'name {cen_atom}')
        atom2_atoms = u.select_atoms(f'name {sur_atom}')
        atom1_indices = atom1_atoms.indices
        atom2_indices = atom2_atoms.indices
        coordination_numbers = []
        neighbor_indices = []
        for atom in atom1_atoms:
            ns = NeighborSearch.AtomNeighborSearch(atom1_atoms)
            atoms_distances = distances.distance_array(atom.position, atom2_atoms.positions)
            distance_cutoff = cutoff
            coordination_number = np.sum(atoms_distances <= distance_cutoff)
            coordination_numbers.append(coordination_number)
            indices = np.where(atoms_distances <= distance_cutoff)[1]
            neighbor_indices.append(list(atom2_atoms.indices[indices]))
        my_bar.progress((index + 1) / len(uploaded_file_paths), text=progress_text)
        df[index+1] = coordination_numbers
        neighbor_df[index+1] = neighbor_indices
    df.insert(0, 'Atom Index', atom1_indices + 1)
    neighbor_df.insert(0, 'Atom Index', atom1_indices + 1)
    time.sleep(1)
    coordination_csv = convert_df(df)
    neighbor_csv = convert_df(neighbor_df)
    st.download_button(
        label="Download coordination as CSV",
        data=coordination_csv,
        file_name=f'{cen_atom}_Cn.csv',
        mime='text/csv',
    )
    st.download_button(
        label="Download neighbor as CSV",
        data=neighbor_csv,
        file_name=f'{cen_atom}_neighbor_{sur_atom}_indices.csv',
        mime='text/csv',
    )
else:
    st.write('Goodbye')