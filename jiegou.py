# -*- coding: utf-8 -*-
# @Time : 2023/4/14 15:53
# @Author : lza
# @File : SCv2.py
# @Software: PyCharm


import copy
import re
import os
import pandas as pd
import time
from collections import Counter
from functools import reduce
import streamlit as st


# 删除原子后的换行符
def del_break(sentence):
    if '\n' in sentence:
        sentence = re.split(r'\n', sentence)[0]
    return sentence


def add_values(key, value, dic):
    elements = ['Si', 'Al', 'Mg', 'Y', 'B']               # 添加Y/Mg
    temp = []
    if re.split('[0-9]+', value)[0] in elements:
        if key not in dic:
            temp.append(value)
            dic[key] = temp
        else:
            if value not in dic[key]:
                dic[key].append(value)


# 查找区分阳离子
def find_polytope(key, value, dic):
    temp = []
    if key not in dic:
        temp.append(value)
        dic[key] = temp
    else:
        if value not in dic[key]:
            dic[key].append(value)
    return dic


# 复制并返回带有所有成键原子的列表，形如[[Si1, O2], [Si1, O3], ...]
# 返回一个氧字典，形如{"O2": ['Si1', 'Al3'], "O4": ['Si1', 'Al3'], ...}
# 返回所有阳离子字典列表，形如[{"Si1": ['O2', 'O3'], ...}, {"Al4": ['O2', 'O3'], ...}, ...]
def find_atoms(filename):
    atoms = []                                  # 存放成键的原子对
    all_oxygen = []                             # 存放所有氧原子
    oxygen_cation = {}                          # 存放与阳离子相连的氧
    elements = ['Si', 'Al', 'Mg', 'Y', 'B']          # 多面体元素，依据实际添加删除
    si_dic = {}
    al_dic = {}
    mg_dic = {}
    b_dic = {}                                  # 增加元素需要同步增加
    ca_dic = {}
    y_dic = {}
    elements_dics_list = [si_dic, al_dic, mg_dic, y_dic, b_dic]    # 存放多面体字典
    file = open(filename)
    lines = file.readlines()
    index_str = "_ccdc_geom_bond_type\n"
    start_index = lines.index(index_str) + 1
    lines = lines[start_index:]
    for line in lines:
        fst = re.split(r' +', line)[0]
        snd = re.split(r' +', line)[1]
        fst = del_break(fst)
        snd = del_break(snd)
        temp = [fst, snd]
        atoms.append(temp)
        if 'O' in fst:
            add_values(fst, snd, oxygen_cation)
            if 'Ca' in snd:
                find_polytope(snd, fst, ca_dic)
            try:
                ele = re.split(r'[0-9]+', snd)[0]
                index = elements.index(ele)
                find_polytope(snd, fst, elements_dics_list[index])
            except ValueError:
                pass
            if fst not in all_oxygen:       # 查找所有氧原子
                all_oxygen.append(fst)
        elif 'O' in snd:
            add_values(snd, fst, oxygen_cation)
            if 'Ca' in fst:
                find_polytope(fst, snd, ca_dic)
            try:
                ele = re.split(r'[0-9]+', fst)[0]
                index = elements.index(ele)
                find_polytope(fst, snd, elements_dics_list[index])
            except ValueError:
                pass
            if snd not in all_oxygen:
                all_oxygen.append(snd)
    file.close()
    free_oxygen = len(all_oxygen) - len(oxygen_cation)     # 游离氧
    df = pd.DataFrame([len(all_oxygen), free_oxygen], index=['Total Oxygen', 'Free Oxygen'], columns=['number'])
    return atoms, oxygen_cation, elements_dics_list, df, ca_dic


def judge_cation(ele_dics):
    num = len(ele_dics)
    si_dic = {}
    al_dic = {}
    mg_dic = {}
    b3_dic = {}
    b4_dic = {}
    y_dic = {}
    poly_ele_dics_list = [si_dic, al_dic, mg_dic, y_dic, b3_dic, b4_dic]        # b3_dic, b4_dic位置不要变动
    for i in range(num):
        ele_dic = ele_dics[i]                         # 依次取出存放Si, Al, Mg, B的字典
        for ele in ele_dic:                           # 这一部分以后可以加入判断双键
            if 'B' not in ele:
                if len(ele_dic[ele]) == 4:
                    poly_ele_dics_list[i][ele] = ele_dic[ele]
            elif 'B' in ele:
                if len(ele_dic[ele]) == 3:
                    poly_ele_dics_list[i][ele] = ele_dic[ele]
                elif len(ele_dic[ele]) == 4:
                    poly_ele_dics_list[i+1][ele] = ele_dic[ele]
    return poly_ele_dics_list


# 计算桥氧非桥氧
# 分别返回桥氧和非桥氧的字典
def devide_bo_nbo(oxy_cat, ele_dics):
    sin_oxygen = []                                      # 一配位氧
    dou_oxygen = []                                      # 二配位氧
    tri_oxygen = []                                      # 三配位氧
    qua_oxygen = []                                      # 四配位氧
    pen_oxygen = []                                      # 五配位氧
    oxy_cn = [sin_oxygen, dou_oxygen, tri_oxygen, qua_oxygen, pen_oxygen]
    oxy_cn_str = ["Sin-Oxygen", "Dou-Oxygen", "Tri-Oxygen", "Qua-Oxygen", "Pen-Oxygen"]
    nbo = {}                                            # 非桥氧字典
    bo = {}                                             # 桥氧字典
    dou_bo = {}
    tri_bo = {}
    for oxy in oxy_cat:
        temp = []
        count = 0
        if len(oxy_cat[oxy]) <= 5:
            oxy_cn[(len(oxy_cat[oxy])-1)].append(oxy)
        for cat in oxy_cat[oxy]:
            if 'B' not in cat:                          # 判断氧连接了几个[Mg04], [Si04], [Al04]
                for ele_dic in ele_dics:
                    if cat in ele_dic:
                        temp.append(cat)
                        count += 1
                        continue
            else:
                if cat in ele_dics[-1]:
                    temp.append(cat)
                    count += 1
                elif cat in ele_dics[-2]:
                    temp.append(cat)
                    count += 1
        if count >= 2:
            bo[oxy] = temp
            if count == 2:
                dou_bo[oxy] = temp
            else:
                tri_bo[oxy] = temp
        else:
            nbo[oxy] = temp
    df = pd.DataFrame([len(bo), len(nbo)], index=['Bridge Oxygen', 'Non Bridge Oxygen'], columns=['number'])
    df2 = pd.DataFrame()
    for i in range(len(oxy_cn)):
        count = [[], [], [], [], [], []]
        count1 = []
        for item in oxy_cn[i]:
            try:
                num = len(bo[item])
                count[num].append(item)
            except KeyError:
                count[0].append(item)
        for item in count:
            count1.append(len(item))
        df1 = pd.DataFrame([count1], columns=['O0', 'O1', 'O2', 'O3', 'O4', 'O5'], index=[oxy_cn_str[i]])
        df2 = pd.concat([df1, df2])
    df2.loc['Total'] = df2.apply(lambda x: x.sum())
    df2['Total'] = df2.apply(lambda x: x.sum(), axis=1)
    return bo, nbo, df, df2, dou_bo, tri_bo


# 判断氧连的多面体类型
def judge_oxy_cation(co_oxy_dic, b3):
    poly_dict = {}
    for oxygen in co_oxy_dic:
        temp = []  # save polyhedra name
        for ion in co_oxy_dic[oxygen]:
            element = re.split('[0-9]+', ion)[0]
            if element == 'B':
                element = 'b3' if ion in b3 else 'b4'
            temp.append(element)
        temp = sorted(temp)
        name = reduce(lambda x, y: x+'-O-'+y, temp)
        if name in poly_dict:
            poly_dict[name].append(oxygen)
        else:
            poly_dict[name] = [oxygen]
    df = pd.DataFrame()
    for _ in poly_dict:
        data = {
            'name': pd.Series(_),
            'number': pd.Series(poly_dict[_].__len__())
        }
        df1 = pd.DataFrame(data=data)
        df = pd.concat([df, df1], ignore_index=True)
    return df


# 判断Qn
def cal_qn(ele_dic, bri_oxy):
    ele_name = re.split(r'[0-9]+', list(ele_dic)[0])[0]
    if ele_name == "B":
        if len(ele_dic[list(ele_dic)[0]]) == 3:
            ele_name = "B3"
        else:
            ele_name = "B4"
    qn_dic = {"Q_oth": [], "Q0": [], "Q1": [], "Q2": [], "Q3": [], "Q4": [], "Q5": []}
    for cat in ele_dic:
        count = 0
        for oxy in ele_dic[cat]:
            if oxy in bri_oxy:
                count += 1
        if count <= 5:
            qn_dic[list(qn_dic)[count+1]].append({cat: ele_dic[cat]})
        else:
            qn_dic["Q_oth"].append({cat: ele_dic[cat]})
    dic = {ele_name: qn_dic}
    dic = copy.deepcopy(dic)
    for qn in qn_dic:
        qn_dic[qn] = len(qn_dic[qn])
    df = pd.DataFrame(qn_dic, index=[ele_name])
    return df, dic


def cal_cn(ele_dic):
    cn_dic = {
        "CN0": [], "CN1": [], "CN2": [], "CN3": [],
        "CN4": [], "CN5": [], "CN6": [], "CN7": [],
        "CN8": [], "CN_other": []
    }
    if len(ele_dic) != 0:
        ele_name = re.split(r'[0-9]+', list(ele_dic)[0])[0]
        for ele in ele_dic:
            index1 = len(ele_dic[ele])
            try:
                cn_dic[list(cn_dic)[index1]].append(1)
            except IndexError:
                cn_dic["CN_other"].append(1)
        for item in cn_dic:
            cn_dic[item] = len(cn_dic[item])
        df = pd.DataFrame(cn_dic, index=[ele_name])
        return df


# 查看修饰体(Mg, Ca, Na, Al)并返回一个字典
def cal_mod_boron(atoms):
    elements = ['Na', 'Ca', 'Mg', 'Al', 'Li']             # 根据需求添加
    poly_element = ['Si', 'B', 'Al', 'Mg']          # 多面体元素列表, 根据需求添加
    na_dic, ca_dic, mg_dic, al_dic, li_dic = {}, {}, {}, {}, {}
    mod_oxy_dic = {}
    ele_mod_dics_list = [na_dic, ca_dic, mg_dic, al_dic, li_dic]
    for atom in atoms:
        fst = atom[0]
        snd = atom[1]
        if re.split(r'[0-9]+', fst)[0] in poly_element:         # 统计形成体元素与修饰体连接情况
            modify_name = re.split(r'[0-9]+', snd)[0]
            try:
                i = elements.index(modify_name)
                find_polytope(snd, fst, ele_mod_dics_list[i])
            except ValueError:
                pass
        elif re.split(r'[0-9]+', snd)[0] in poly_element:
            modify_name = re.split(r'[0-9]+', fst)[0]
            try:
                i = elements.index(modify_name)
                find_polytope(fst, snd, ele_mod_dics_list[i])
            except ValueError:
                pass
        if 'O' in fst:
            modify_name = re.split(r'[0-9]+', snd)[0]
            if modify_name in elements:
                find_polytope(fst, snd, mod_oxy_dic)
        elif 'O' in snd:
            modify_name = re.split(r'[0-9]+', fst)[0]
            if modify_name in elements:
                find_polytope(snd, fst, mod_oxy_dic)
    return ele_mod_dics_list, mod_oxy_dic


# 排除成四面体的Mg/Al
def del_poly_mod_ele(poly_ele_dics, ele_mod_dics):
    index1, index2 = None, None
    index11, index22 = None, None
    al_dic = {}
    mg_dic = {}
    for i in range(len(poly_ele_dics)):     # 找到Al, Mg 在多面体字典列表里的索引
        poly_ele_list = list(poly_ele_dics[i])
        if len(poly_ele_list) != 0:
            if "Al" in poly_ele_list[0]:
                index1 = i              # Al 在多面体字典列表里的序列
                continue
            elif "Mg" in poly_ele_list[0]:
                index2 = i              # Mg 在多面体字典列表的序列
                continue
    for i in range(len(ele_mod_dics)):      # 找到Al, Mg 在编辑体字典列表里的索引
        ele_mod_list = list(ele_mod_dics[i])
        if len(ele_mod_list) != 0:
            if "Al" in ele_mod_list[0]:
                index11 = i             # Al 在编辑体字典列表里的序列
                continue
            elif "Mg" in ele_mod_list[0]:
                index22 = i             # Mg 在编辑体字典列表里的序列
                continue
    if index1 and index11 is not None:
        for item in ele_mod_dics[index11]:          # 有些修饰体没有与多面体成键
            if item in poly_ele_dics[index1]:
                continue
            else:
                al_dic[item] = ele_mod_dics[index11][item]
    if index2 and index22 is not None:
        for item in ele_mod_dics[index22]:
            if item in poly_ele_dics[index2]:
                continue
            else:
                mg_dic[item] = ele_mod_dics[index22][item]
    if len(al_dic) != 0:
        ele_mod_dics[index11] = al_dic
    if len(mg_dic) != 0:
        ele_mod_dics[index22] = mg_dic
    return ele_mod_dics


# 分析多面体与修饰体的连接情况
def ana_poly_mod(ele_mod_dics, poly_ele_dics):
    poly_mod_dic = {}
    poly_name_list = []
    for item in poly_ele_dics:
        if len(item) != 0:
            poly_name = re.split(r"[0-9]+", list(item)[0])[0]
            if 'B' in poly_name and "B3" not in poly_name_list:   # 依托于前面的多面体顺序，B3/B4顺序不能变，即先B3再B4
                poly_name = "B3"
            elif 'B' in poly_name and 'B3' in poly_name_list:
                poly_name = "B4"
            poly_name_list.append(poly_name)
            continue
    for ele_mod_dic in ele_mod_dics:
        if len(ele_mod_dic) != 0:
            mod_ele_name = re.split(r"[0-9]+", list(ele_mod_dic)[0])[0]     # 编辑体名字
            for item in ele_mod_dic:
                for ele in ele_mod_dic[item]:
                    ele_name = re.split(r"[0-9]+", ele)[0]                  # 相连的多面体名字
                    if 'B' in ele_name:
                        index1 = poly_name_list.index("B3")
                        index2 = poly_name_list.index("B4")
                        if ele in poly_ele_dics[index1]:
                            ele_name = "B3"
                            key = ele_name + '-' + mod_ele_name
                            if key not in poly_mod_dic:
                                poly_mod_dic[key] = [ele]
                            else:
                                poly_mod_dic[key].append(ele)
                        elif ele in poly_ele_dics[index2]:
                            ele_name = "B4"
                            key = ele_name + '-' + mod_ele_name
                            if key not in poly_mod_dic:
                                poly_mod_dic[key] = [ele]
                            else:
                                poly_mod_dic[key].append(ele)
                    else:
                        try:
                            index1 = poly_name_list.index(ele_name)
                            if ele in poly_ele_dics[index1]:
                                key = ele_name + '-' + mod_ele_name
                                if key not in poly_mod_dic:
                                    poly_mod_dic[key] = [ele]
                                else:
                                    poly_mod_dic[key].append(ele)
                        except ValueError:
                            continue
            continue
    df_poly_mod = pd.DataFrame()
    for item in poly_mod_dic:
        data = {
            'name': pd.Series(item),
            'number': pd.Series(poly_mod_dic[item].__len__())
        }
        df1 = pd.DataFrame(data=data)
        df_poly_mod = pd.concat([df_poly_mod, df1], ignore_index=True)
    return df_poly_mod


# 删除由cal_mod_boron(atoms)函数返回的 mod_oxy_dics_list内的修饰体
def del_mod_oxy_poly(modify_oxygen_dic, polytopes_dics_list):
    mod_oxy_dic = {}
    poly_name = []
    mod_ele = ['Na', 'Ca', 'Mg', 'Li']
    for item in polytopes_dics_list:
        try:
            ele_name = re.split(r'[0-9]+', list(item)[0])[0]
            poly_name.append(ele_name)
        except IndexError:
            continue
    for item in modify_oxygen_dic:
        temp = []
        for ele in modify_oxygen_dic[item]:
            ele_name = re.split(r'[0-9]+', ele)[0]
            try:
                if ele not in polytopes_dics_list[poly_name.index(ele_name)]:
                    temp.append(ele)
            except ValueError:
                if ele_name in mod_ele:
                    temp.append(ele)
        modify_oxygen_dic[item] = temp
    for item in modify_oxygen_dic:
        if len(modify_oxygen_dic[item]) == 0:
            continue
        else:
            mod_oxy_dic[item] = modify_oxygen_dic[item]
    return mod_oxy_dic


# 统计形成体内的氧与修饰体相连的情况
def ana_mod_poly_oxy(mod_oxy_dic, qn_list):
    df = pd.DataFrame()
    for polyhedron in qn_list:
        for item in polyhedron:
            ele_name = item
            for qn in polyhedron[item]:
                if len(polyhedron[item][qn]):
                    bridge_oxygen_number = qn
                    temp_dic = {}
                    for polytope_oxygen in polyhedron[item][qn]:
                        modify_temp = []
                        # 下面遍历是找一个形成体内氧于编辑体的连接情况
                        for oxygens in polytope_oxygen.values():
                            modify_duibi_list = []  # 存放连接的修饰体
                            for oxygen in oxygens:
                                try:
                                    for modify in mod_oxy_dic[oxygen]:
                                        if modify not in modify_duibi_list:
                                            modify_duibi_list.append(modify)
                                            modify_name = re.split(r"[0-9]+", modify)[0]
                                            modify_temp.append(modify_name)
                                except KeyError:
                                    pass

                        if modify_temp.__len__():
                            modify_dic = dict(Counter(modify_temp))
                            temp = []
                            for ion in modify_dic:
                                name = str(modify_dic[ion]) + ion
                                temp.append(name)
                            key = ""
                            for i in sorted(temp):
                                key = key + '-' + i
                            if key not in temp_dic:
                                temp_dic[key] = [polytope_oxygen]
                            else:
                                temp_dic[key].append(polytope_oxygen)

                    for i in temp_dic:
                        temp_dic[i] = len(temp_dic[i])
                    df1 = pd.DataFrame(data=temp_dic, index=[ele_name+bridge_oxygen_number])
                    df = pd.concat([df, df1])

    df_oxy_mod = df.fillna(0)
    return df_oxy_mod


def cal_ave_mod(df_oxy_mod, df_qn):

    cols = [i for i in df_oxy_mod.columns if i != 'Total']
    number = []

    for col in cols:
        elements = re.split(r'-', col)
        elements = [i for i in elements if i != '']
        temp = []
        for element in elements:
            num = re.split(r'[a-zA-Z]+', element)[0]
            ele = re.split(r'[0-9]+', element)[-1]
            temp.append((ele, num))
        number.append(temp)

    dic = {}
    for ind in df_oxy_mod.index:
        dic.setdefault(ind, {})
        temp = [i for i in df_oxy_mod.loc[ind]]
        for i in range(temp.__len__()):
            ele = number[i]
            if temp[i] == 0:
                continue
            else:
                for t in ele:
                    key = t[0]
                    multi = t[1]
                    dic[ind].setdefault(key, 0)
                    num = int(multi)*int(temp[i])
                    before_sum = dic[ind][key]
                    after_sum = before_sum + num
                    dic[ind][key] = after_sum

    df_cal_ave = (pd.DataFrame(data=dic)).T
    df_cal_ave = df_cal_ave.fillna(0)

    temp = {}
    for name in df_cal_ave.index:
        ele = name[:2]
        qn = name[2:]
        temp[name] = df_qn.loc[ele][qn]

    df_oxy_mod['Total'] = df_oxy_mod.apply(lambda x: x.sum(), axis=1)

    temp1, temp2 = [], []           # temp1 save num/Qn, temp2 save num/Total
    for ind in df_cal_ave.index:
        temp1.append([df_cal_ave.loc[ind][col]/int(temp[ind]) for col in df_cal_ave.columns])
        temp2.append([df_cal_ave.loc[ind][col]/int(df_oxy_mod.loc[ind]['Total']) for col in df_cal_ave.columns])

    name1 = list(df_cal_ave.columns)
    name2 = [i+'_ave_total' for i in name1]
    name1 = [i + '_ave_qn' for i in name1]
    df_temp = pd.DataFrame(data=temp1, columns=name1, index=df_cal_ave.index)
    df_temp1 = pd.DataFrame(data=temp2, columns=name2, index=df_cal_ave.index)
    df_cal_ave = pd.concat([df_cal_ave, df_temp, df_temp1], axis=1)

    return df_cal_ave


def link_single_modify(df_oxy_mod):

    cols = list(df_oxy_mod.columns)
    cols = [i for i in cols if len(re.findall('-', i)) == 1]

    df_single = pd.DataFrame(index=df_oxy_mod.index, columns=cols)
    df_single = df_single.sort_index(axis=1)

    temp_col = []
    for col in cols:
        df_single[col] = df_oxy_mod[col]
        num = re.split(r'[a-zA-Z]+', col)[0]
        num = int(re.split(r'-', num)[1])
        name = re.split(r'[0-9]+', col)[-1]
        temp_col.append((name, num))

    temp = []
    temp1 = []
    mod_name = []
    former_temp = []
    for ind in df_single.index:
        dic, dic1 = {}, {}
        dic.setdefault(ind, {})
        dic1.setdefault(ind, {})
        temp_index = [i for i in df_single.loc[ind]]
        for i in range(len(temp_index)):
            multi = temp_col[i][1]
            name = temp_col[i][0]
            mod_name.append(name if name not in mod_name else None)
            dic[ind].setdefault(name, 0)
            dic1[ind].setdefault(name, 0)
            num = int(multi*int(temp_index[i]))
            dic[ind][name] = dic[ind][name] + num
            dic1[ind][name] = dic1[ind][name] + int(temp_index[i])
        for item in dic:
            temp.append([i for i in list(dic[item].values())])
        for item in dic1:
            temp1.append([i for i in list(dic1[item].values())])
    mod_name = [i for i in mod_name if i is not None]

    # df_single['Total'] = df_single.apply(lambda x: x.sum(), axis=1)

    df_temp = pd.DataFrame(data=temp, columns=mod_name, index=df_single.index)
    df_single = pd.concat([df_single, df_temp], axis=1)

    name2 = [i + '_tol' for i in mod_name]
    df_temp = pd.DataFrame(data=temp1, columns=name2, index=df_single.index)
    df_single = pd.concat([df_single, df_temp], axis=1)

    name1 = [i+'_ave' for i in mod_name]

    for i in range(len(mod_name)):
        value = df_single[mod_name[i]]/df_single[name2[i]]
        df_single.insert(loc=len(df_single.columns), column=name1[i], value=value)

    for item in name2:
        df_single.pop(item)

    df_single = df_single.fillna(0)

    return df_single


def main(filename):
    start_time = time.time()

    file = re.split(r'\.cif', filename)[0]
    file = file + '.xlsx'
    print(f"文件保存名称:{file}")

    atom_pairs, oxy_cat, origin_dics_list, df_oxygen, ca_dic = find_atoms(filename)
    poly_dics_list = judge_cation(origin_dics_list)

    poly_dics_list[-3] = {}  # 需要Y注释这里
    poly_dics_list[-4] = {}  # 需要Mg注释这里

    bri_oxy, non_bri_oxy, df_temp, df_polyconnected, dbo, tbo = devide_bo_nbo(oxy_cat, poly_dics_list)
    print(non_bri_oxy)
    df_oxygen = pd.concat([df_temp, df_oxygen])
    with pd.ExcelWriter(file) as writer:
        df_oxygen.to_excel(writer, sheet_name='氧情况')
        df_polyconnected.to_excel(writer, sheet_name="氧配位以及连接多面体情况")
    print('success save oxygen CN!')

    boron3 = []
    for ele in poly_dics_list[-2]:
        boron3.append(ele)
    dbo_df_cations = judge_oxy_cation(dbo, boron3)
    tbo_df_cations = judge_oxy_cation(tbo, boron3)

    df_qn = pd.DataFrame()

    qn_list = []
    for ele_dic in poly_dics_list:
        if len(ele_dic) != 0:
            df_temp, dic = cal_qn(ele_dic, bri_oxy)
            df_qn = pd.concat([df_qn, df_temp])
            qn_list.append(dic)
    df_qn['Total'] = df_qn.apply(lambda x: x.sum(), axis=1)
    df_qn.loc['Total'] = df_qn.apply(lambda x: x.sum())
    index1 = df_qn.axes[0]
    nbo_t = []
    for item in index1:
        nbo_t1 = 0
        for i in range(1, 6):
            nbo_t1 += df_qn.loc[item][i]*(5-i)
        nbo_t1 /= df_qn.loc[item][-1]
        nbo_t.append(nbo_t1)
    df_temp = pd.DataFrame(nbo_t, index=index1, columns=['NBO/T'])
    df_qn = pd.concat([df_qn, df_temp], axis=1)

    df_cn = pd.DataFrame()
    for ele_dic in origin_dics_list:
        if len(ele_dic) != 0:
            df1 = cal_cn(ele_dic)
            df_cn = pd.concat([df_cn, df1])
    df_temp = cal_cn(ca_dic)
    df_cn = pd.concat([df_cn, df_temp])
    df_cn['Total'] = df_cn.apply(lambda x: x.sum(), axis=1)
    df_cn.loc['Total'] = df_cn.apply(lambda x: x.sum())

    modify_ele_dics_list, modify_oxygen_dics_list = cal_mod_boron(atom_pairs)
    modify_ele_dics_list = del_poly_mod_ele(poly_dics_list, modify_ele_dics_list)
    df_poly_mod = ana_poly_mod(modify_ele_dics_list, poly_dics_list)

    modify_oxygen_dic = del_mod_oxy_poly(modify_oxygen_dics_list, poly_dics_list)
    df_oxy_mod = ana_mod_poly_oxy(modify_oxygen_dic, qn_list)
    df_oxy_mod = df_oxy_mod.sort_index(axis=1)
    df_ave_mod = cal_ave_mod(df_oxy_mod, df_qn)

    df_single_mod = link_single_modify(df_oxy_mod)


    with pd.ExcelWriter(file, mode="a") as writer:
        dbo_df_cations.to_excel(writer, sheet_name="多面体之间连接情况(DBO)")
        tbo_df_cations.to_excel(writer, sheet_name="多面体之间连接情况(TBO)")
        print("success save polytopes connected!")

        df_qn.to_excel(writer, sheet_name="Qn")
        print("success save Qn!")

        df_cn.to_excel(writer, sheet_name='CN')
        print('success save CN!')

        df_poly_mod.to_excel(writer, sheet_name="修饰体直连形成体离子情况")
        print('success save Number Of Modifiers Connect With Polytopes, straightly!')

        df_oxy_mod.to_excel(writer, sheet_name="修饰体与连接体的氧连接情况")
        df_ave_mod.to_excel(writer, sheet_name="修饰体统计0")
        df_single_mod.to_excel(writer, sheet_name="修饰体统计1")
        print("success save!")
    print("用时: %.4f s" % (time.time() - start_time))


if __name__ == "__main__":
    title = st.text_input('文件地址', '请输入')
    st.write('文件地址 is', title)
    if st.button('统计配位信息'):
        main(rf'{title}')
    else:
        st.write('xlsx文件将会保存在文件地址中')



