import pandas as pd
import numpy as np
import os
from itertools import combinations, combinations_with_replacement, product
import string
import matplotlib.pyplot as plt


def fix_time_label(data):
    # TODO
    return complete_data_label(data)

def check_data_cols(data):
    
    cols = data.columns.values
    index = data.index.values
    
    print(np.all(cols == index))

def complete_data_label(data):

    cols = data.columns.values
    complete_labels = [complete_label(item) for item in cols]
    complete_labels

    data.index = complete_labels
    data.columns = complete_labels
    
    return data

def get_same_row_col_labels(data):
    
    cols = set(data.columns)
    rows = set(data.index)
    row_col_labels = sorted(list(cols.intersection(rows)))
    
    return data.loc[row_col_labels, row_col_labels]

def get_cell_diagonal_labels(n):
    
    numbers = np.arange(2, n+2)
    letters = list(string.ascii_uppercase)[1:n+1]
    
    labels = [letters[i] + str(numbers[i]) for i in range(n)]
    
    return labels

def complete_label(label):
    
    label_comps = label.split("_")
    
    if len(label_comps) == 3:
        label_comps.insert(2, '0h')

    new_label = '_'.join(label_comps)
    
    return new_label

def get_var(label, idx):
    
    label_comps = label.split('_')
    
    return label_comps[idx]

def get_var_row_subset_df(var_name, var_idx, data):
    
    is_var_row = [(var_name == get_var(name, var_idx)) for name in data.index]
    
    return data.loc[is_var_row, :].copy()

def get_var_pair_df(var_pair, data, exclude=False):
    
    
    cols = data.columns
    rows = data.index
    
    cols_comps = [item.split('_') for item in cols]
    
    if exclude:
        is_var_2 = [(var_pair[1] not in item.split('_'))  for item in cols]
        is_var_1 = [(var_pair[0] not in item.split('_'))  for item in rows]
    else:
        is_var_2 = [(var_pair[1] in item.split('_'))  for item in cols]
        is_var_1 = [(var_pair[0] in item.split('_'))  for item in rows]
    
    
    var_pair_subset = data.loc[is_var_1, is_var_2]    
    
    return var_pair_subset

def get_var_list(var_idx, data):
    
    var_list = sorted(set([get_var(name, var_idx) for name in data.columns]))
    
    return var_list

def get_var_pair(var_idx, data, include_permutation=False):
    
    var_list = sorted(set([get_var(name, var_idx) for name in data.columns]))

    if include_permutation:
        var_pairs = list(product(var_list, repeat=2))
    else:
        var_pairs = list(combinations(var_list, 2))

        for item in var_list:
            self_self = item, item
            var_pairs.append(self_self)
    
    return var_pairs
    


def get_var_commons(ref_var_idx, common_var_idx, data ):
    
    ref_var_list = get_var_list(ref_var_idx, data)
    
    ref_vars_vars = {}
    for ref_var in ref_var_list:

        ref_var_vars = set([get_var(col_name, common_var_idx)  for col_name in data.columns if ref_var in col_name])
        ref_vars_vars[ref_var] = ref_var_vars
    
    common_vars = list(set.intersection(*list(ref_vars_vars.values())))
    
    return common_vars


def filter_var_subset_df(var_list, var_idx, data):
    
    cols = data.columns
    
    has_var_list_label =  [get_var(name, var_idx) in var_list for name in cols]
    
    df_var_list_subset = data.loc[has_var_list_label, has_var_list_label]
    
    return df_var_list_subset
    

def get_likelihood(x, mu=0, sigma=1):
    var = sigma **2
    L = (1 / (np.sqrt(2 * np.pi * var))) * \
    np.exp(-(x - mu)**2 / (2*var))
    
    return L

def geo_mean(iterable):
    a = np.array(iterable)
    return (a**(1.0/len(a))).prod()

def get_strain_var_likelihood(var_str_idx, var_2_idx, data,  var2_pair=None, geom_normalize=True, use_common_var_2=False, debug=False):
    
    df = data.copy()
    
    # filter df so as to have only organs common across all strains
    if use_common_var_2 is True:
        common_var_2s = get_var_commons(ref_var_idx=var_str_idx, common_var_idx=var_2_idx, data=df)
        df = filter_var_subset_df(var_list=common_var_2s, var_idx=var_2_idx, data=df).copy()
    
    # check if there is a variable 2 to filter 
    if var2_pair is not None:
        df = get_var_pair_df(var2_pair, df)

    var_1_pairs = get_var_pair(var_idx=var_str_idx, data=df, include_permutation=False)
#     print(var_1_pairs)
    var_2_pairs_combination = get_var_pair(var_idx=var_2_idx, data=df, include_permutation=False)
    var_2_pairs_permuted = get_var_pair(var_idx=var_2_idx, data=df, include_permutation=True)
#     print(var_2_pairs)
    
    var_1_list = get_var_list(var_idx=var_str_idx, data=df)

    var_var_df = pd.DataFrame(index=var_1_list, columns=var_1_list)
    
    if debug:
        var_pairs_df = pd.DataFrame(index=var_2_pairs_permuted, columns=var_1_pairs)
        var_pairs_df.columns = var_pairs_df.columns.map('-'.join)
        var_pairs_df.index = var_pairs_df.index.map('-'.join)

        var_pairs_mean_df =  pd.DataFrame(index=var_pairs_df.index, columns=var_pairs_df.columns)



    for i, var_1_pair in enumerate(var_1_pairs):
        
        corr_df_var_1_x = get_var_pair_df((var_1_pair[0], var_1_pair[0]), data=df)
        corr_df_var_1_y = get_var_pair_df((var_1_pair[1], var_1_pair[1]), data=df)


        if var_1_pair[0] != var_1_pair[1]:
            var_2_pairs = var_2_pairs_permuted
        else:
            var_2_pairs = var_2_pairs_combination
        
        var_2_lhs = []
        for j, var_2_pair in enumerate(var_2_pairs):
            
            corr_df_var_1_x_var_2 = get_var_pair_df((var_2_pair[0], var_2_pair[1]), data=corr_df_var_1_x).values
            corr_df_var_1_y_var_2 = get_var_pair_df((var_2_pair[0], var_2_pair[1]), data=corr_df_var_1_y).values
            

            # take upper diagonal if same var_2 (e.g same organs)
            if var_2_pair[0] == var_2_pair[1]:
                corr_df_var_1_x_var_2 = corr_df_var_1_x_var_2[np.tril_indices_from(corr_df_var_1_x_var_2)]
                corr_df_var_1_y_var_2 = corr_df_var_1_y_var_2[np.tril_indices_from(corr_df_var_1_y_var_2)]
            

            if (corr_df_var_1_x_var_2.size >= 2) and (corr_df_var_1_y_var_2.size >= 2):

                # row is reference, row is mean
                mu = corr_df_var_1_x_var_2.mean()
                sigma = corr_df_var_1_x_var_2.std(ddof=1)  # sigma can be zero if all values are same
                x = corr_df_var_1_y_var_2
                
                var_1_pair_to_check = (('AJ', 'BL'), ('BL', 'BL')); var_2_pair_to_check = (('bone', 'brain'), ('bone', 'bone'))
                
                if var_1_pair in var_1_pair_to_check and var_2_pair in var_2_pair_to_check:
                    print (f"mu for {var_1_pair} and {var_2_pair} = {mu} and std = {sigma}")
                    
                if sigma != 0:
                    var_2_L = get_likelihood(x, mu, sigma)
                    
                
                    if geom_normalize == True:
                # find the n_root to keep likelihood computation uniform for different n_L
                        var_2_L = var_2_L ** (1/var_2_L.size)

                    var_2_L = var_2_L.prod()
                    var_2_lhs.append(var_2_L)

                    if debug:
                        var_pairs_df.loc['-'.join(var_2_pair), '-'.join(var_1_pair)] = var_2_L
                        var_pairs_mean_df.loc['-'.join(var_2_pair), '-'.join(var_1_pair)] = (mu, sigma)

    
    
            
        var_2_lhs = np.array(var_2_lhs)
        var_2_lhs = var_2_lhs[~np.isnan(var_2_lhs)]

        lh =  np.mean(var_2_lhs)
        var_var_df.loc[var_1_pair[0], var_1_pair[1]] = lh
    if debug:
        print('== lh values== \n',var_pairs_df)
        print('== mean and std. values== \n',var_pairs_mean_df)  
    return var_var_df

def get_organ_var_likelihood(var_str_idx, var_2_idx, data,  var2_pair=None, geom_normalize=True, use_common_var_2=False, debug=False):
    
    df = data.copy()
    
    # filter df so as to have only organs common across all strains
    if use_common_var_2 is True:
        common_var_2s = get_var_commons(ref_var_idx=var_2_idx, common_var_idx=var_str_idx, data=df)
        df = filter_var_subset_df(var_list=common_var_2s, var_idx=var_str_idx, data=df).copy()
        
        if debug:
            print(common_var_2s)

    # check if there is a variable 2 to filter 
    if var2_pair is not None:
        df = get_var_pair_df(var2_pair, df)

    var_1_pairs = get_var_pair(var_idx=var_str_idx, data=df, include_permutation=False)
    var_2_pairs_combination = get_var_pair(var_idx=var_2_idx, data=df, include_permutation=False)
    var_2_pairs_permuted = get_var_pair(var_idx=var_2_idx, data=df, include_permutation=True)
    
    var_1_list = get_var_list(var_idx=var_str_idx, data=df)

    var_var_df = pd.DataFrame(index=var_1_list, columns=var_1_list)
    
    if debug:
        var_pairs_df = pd.DataFrame(index=var_2_pairs_permuted, columns=var_1_pairs)
        var_pairs_df.columns = var_pairs_df.columns.map('-'.join)
        var_pairs_df.index = var_pairs_df.index.map('-'.join)
        
        var_pairs_mean_df =  pd.DataFrame(index=var_pairs_df.index, columns=var_pairs_df.columns)



    for i, var_1_pair in enumerate(var_1_pairs):
        
        corr_df_var_1_xy = get_var_pair_df((var_1_pair[0], var_1_pair[1]), data=df)
     
        if var_1_pair[0] != var_1_pair[1]:
            var_2_pairs = var_2_pairs_permuted
        else:
          var_2_pairs = var_2_pairs_combination
        # var_2_pairs = var_2_pairs_combination
        
        var_2_lhs = []
        for j, var_2_pair in enumerate(var_2_pairs):

            corr_df_var_1_xy_var_2_x = get_var_pair_df((var_2_pair[0], var_2_pair[0]), data=corr_df_var_1_xy).values
            corr_df_var_1_xy_var_2_y = get_var_pair_df((var_2_pair[1], var_2_pair[1]), data=corr_df_var_1_xy).values
            

             # take upper diagonal if same organs
            if var_1_pair[0] == var_1_pair[1]:
                corr_df_var_1_xy_var_2_x = corr_df_var_1_xy_var_2_x[np.tril_indices_from(corr_df_var_1_xy_var_2_x)]
                corr_df_var_1_xy_var_2_y = corr_df_var_1_xy_var_2_y[np.tril_indices_from(corr_df_var_1_xy_var_2_y)]

            if (corr_df_var_1_xy_var_2_x.size >= 2) and (corr_df_var_1_xy_var_2_y.size >= 2):


                # row is reference, row is mean
                mu = corr_df_var_1_xy_var_2_x.mean()
                sigma = corr_df_var_1_xy_var_2_x.std(ddof=1)  # sigma can be zero if all values are same
                x = corr_df_var_1_xy_var_2_y

                if debug:
                    var_1_pair_to_check = (('lung', 'kidney'), ('kidney', 'kidney')); var_2_pair_to_check = (('AJ', 'BL'), ('BL', 'BL'))
                    
                    if var_1_pair in var_1_pair_to_check and var_2_pair in var_2_pair_to_check:
                        print (f"mu for {var_1_pair} and {var_2_pair} = {mu} and std = {sigma}")
                
                if sigma != 0:
                
                    var_2_L = get_likelihood(x, mu, sigma)
                    if geom_normalize == True:
                # find the n_root to keep likelihood computation uniform for different n_L
                        var_2_L = var_2_L ** (1/var_2_L.size)

                    var_2_L = var_2_L.prod()
                    var_2_lhs.append(var_2_L)

                    if debug:
                        var_pairs_df.loc['-'.join(var_2_pair), '-'.join(var_1_pair)] = var_2_L
                        var_pairs_mean_df.loc['-'.join(var_2_pair), '-'.join(var_1_pair)] = (mu, sigma)

    
            
        var_2_lhs = np.array(var_2_lhs)
        var_2_lhs = var_2_lhs[~np.isnan(var_2_lhs)]
        
        lh =  np.mean(var_2_lhs)
        var_var_df.loc[var_1_pair[0], var_1_pair[1]] = lh
    
    if debug:
        print('== lh values== \n',var_pairs_df)
        print('== mean and std. values== \n',var_pairs_mean_df)
        
    return var_var_df

def get_ranking(raw_df_path, save_path):
    
    data_hr_pair_mat = pd.ExcelFile(raw_df_path)

    with pd.ExcelWriter(os.path.join(save_path), engine='xlsxwriter') as writer:

        workbook = writer.book
        fail_format = workbook.add_format({'bg_color': '#FFC7CE',
                       'font_color': '#9C0006'})
        for i, sheet in enumerate(data_hr_pair_mat.sheet_names):
            sheet_df = data_hr_pair_mat.parse(sheet_name=sheet, index_col=0)
            #sheet_df = sheet_df.iloc[:, :-2]
            rank_df = sheet_df.rank(ascending=False, method='max', axis=0)
            #rank_df['similarity_rank'] = rank_df.mean(axis=1)
            rank_df.to_excel(writer, sheet_name=sheet)

            nrows = rank_df.shape[0]
            ncols = rank_df.shape[1]

            cell_labels = get_cell_diagonal_labels(nrows)

            worksheet = writer.sheets[sheet]
            for j, cell_label in enumerate(cell_labels):



                worksheet.conditional_format(cell_label, {'type':     'cell',
                                           'criteria': 'not equal to',
                                           'value':    1,
                                       'format':   fail_format})
    data_hr_pair_mat.close()
    


    
    df = data.copy()
    
    # filter df so as to have only organs common across all strains
    if use_common_var_2 is True:
        common_var_2s = get_var_commons(ref_var_idx=var_2_idx, common_var_idx=var_str_idx, data=df)
        df = filter_var_subset_df(var_list=common_var_2s, var_idx=var_str_idx, data=df).copy()

    # check if there is a variable 2 to filter 
    if var2_pair is not None:
        df = get_var_pair_df(var2_pair, df)

    var_1_pairs = get_var_pair(var_idx=var_str_idx, data=df, include_permutation=False)
#     var_2_pairs_combination = get_var_pair(var_idx=var_2_idx, data=df, include_permutation=True)
    var_2_pairs_permuted = get_var_pair(var_idx=var_2_idx, data=df, include_permutation=True)
    
    var_1_list = get_var_list(var_idx=var_str_idx, data=df)
#     var_2_list = get_var_list(var_idx=var_2_idx, data=df)

    # dataframe of non-aggregated likelihood
    var_var_df = pd.DataFrame(index=var_1_list, columns=var_1_list)
    
    var_pairs_df = pd.DataFrame(index=var_2_pairs_permuted, columns=var_1_pairs)
    var_pairs_df.columns = var_pairs_df.columns.map('-'.join)
    var_pairs_df.index = var_pairs_df.index.map('-'.join)
    
    var_pairs_mean_df =  pd.DataFrame(index=var_pairs_df.index, columns=var_pairs_df.columns)

    for i, var_1_pair in enumerate(var_1_pairs):
        
        corr_df_var_1_xy = get_var_pair_df((var_1_pair[0], var_1_pair[1]), data=df)

        var_2_pairs = var_2_pairs_permuted
    
        var_2_lhs = []
        for j, var_2_pair in enumerate(var_2_pairs):
            
            corr_df_var_1_xy_var_2_x = get_var_pair_df((var_2_pair[0], var_2_pair[0]), data=corr_df_var_1_xy).values
            corr_df_var_1_xy_var_2_y = get_var_pair_df((var_2_pair[1], var_2_pair[1]), data=corr_df_var_1_xy).values
            

             # take upper diagonal if same organs
            if var_1_pair[0] == var_1_pair[1]:
                corr_df_var_1_xy_var_2_x = corr_df_var_1_xy_var_2_x[np.tril_indices_from(corr_df_var_1_xy_var_2_x)]
                corr_df_var_1_xy_var_2_y = corr_df_var_1_xy_var_2_y[np.tril_indices_from(corr_df_var_1_xy_var_2_y)]

            if (corr_df_var_1_xy_var_2_x.size >= 2) and (corr_df_var_1_xy_var_2_y.size >= 2):


                # row is reference, row is mean
                mu = corr_df_var_1_xy_var_2_x.mean()
                sigma = corr_df_var_1_xy_var_2_x.std(ddof=1)  # sigma can be zero if all values are same
                x = corr_df_var_1_xy_var_2_y

                var_1_pair_to_check = (('lung', 'kidney'), ('kidney', 'kidney')); var_2_pair_to_check = (('AJ', 'BL'), ('BL', 'BL'))
                
                if var_1_pair in var_1_pair_to_check and var_2_pair in var_2_pair_to_check:
                    print (f"mu for {var_1_pair} and {var_2_pair} = {mu} and std = {sigma}")
                
                if sigma != 0:
                
                    var_2_L = get_likelihood(x, mu, sigma)
                    if geom_normalize == True:
                # find the n_root to keep likelihood computation uniform for different n_L
                        var_2_L = var_2_L ** (1/var_2_L.size)
                    
                    var_2_L = var_2_L.prod()
                    var_2_lhs.append(var_2_L)
                    
                    
                    var_pairs_df.loc['-'.join(var_2_pair), '-'.join(var_1_pair)] = var_2_L
                    var_pairs_mean_df.loc['-'.join(var_2_pair), '-'.join(var_1_pair)] = (mu, sigma)

                    
                    
    
            
        var_2_lhs = np.array(var_2_lhs)
        var_2_lhs = var_2_lhs[~np.isnan(var_2_lhs)]
        
        lh = np.mean(var_2_lhs)
        if var_1_pair in [('bone', 'brain')]:
            print(var_2_lhs)
            print(lh)
        var_var_df.loc[var_1_pair[0], var_1_pair[1]] = lh
        
    print('== lh values== \n',var_pairs_df)
    print('== mean and std. values== \n',var_pairs_mean_df)
    return var_var_df