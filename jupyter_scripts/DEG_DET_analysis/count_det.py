#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

def count_rows(filepath):
    """统计文件的行数（排除首行）"""
    with open(filepath, 'r', encoding='utf-8') as f:
        # 跳过表头
        next(f, None)
        return sum(1 for _ in f)

def main():
    root_dir = os.getcwd()
    results = []

    # 遍历当前目录下的每个物种文件夹
    for species in os.listdir(root_dir):
        species_path = os.path.join(root_dir, species)
        if not os.path.isdir(species_path):
            continue

        # 在每个物种文件夹中查找符合命名规则的文件
        for fname in os.listdir(species_path):
            # 只处理两个类型的文件
            if fname.endswith('_DEG.tsv'):
                compare = fname[:-len('_DEG.tsv')]  # 提取 Compare 名称
                det_path = os.path.join(species_path, fname)
                det_count = count_rows(det_path)
                
                # 寻找对应的 nonDEG 文件
                nondeg_fname = compare + '_DTx_nonDEG.tsv'
                nondeg_path = os.path.join(species_path, nondeg_fname)
                if os.path.exists(nondeg_path):
                    nondeg_count = count_rows(nondeg_path)
                else:
                    nondeg_count = 0

                results.append({
                    'species': species,
                    'Compare': compare,
                    'DET': det_count,
                    'DET_non_DEG': nondeg_count
                })

    # 输出结果
    header = ['species', 'Compare', 'DET', 'DET_non_DEG']
    print('\t'.join(header))
    for row in results:
        print(f"{row['species']}\t{row['Compare']}\t{row['DET']}\t{row['DET_non_DEG']}")

if __name__ == '__main__':
    main()

