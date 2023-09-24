def cat_purification(cat):
    cat_complete = cat[(cat['core_score'] != 0) & (cat['notes'] != 'HII')]
    len(cat) - len(cat_complete)
    print('Removed '+str(len(cat) - len(cat_complete))+' sources for the complete catalog')
    cat_medium = cat_complete[cat_complete['notes'] != 'HII?']
    cat_medium = cat_medium[cat_medium['notes'] != 'HII_new?']
    print('Removed '+str(len(cat_complete) - len(cat_medium))+' sources for the medium catalog')
    cat_robust = cat_medium[cat_medium['core_score'] == 2]
    print('Removed '+str(len(cat_medium) - len(cat_robust))+' sources for the robust catalog')
    return cat_complete, cat_medium, cat_robust

def remove_HII_cores(catB3_HII_array, catB6):
    list_of_rows_to_remove = []
    for i in range(len(catB6)):
        B3_name = catB6['B3_match'][i]
        for catB3 in catB3_HII_array:
            if B3_name in catB3['_name']:
                #print('Band 3: '+str(B3_name)+', Band 6: '+str(catB6['_name'][i]))
                list_of_rows_to_remove += [i]
    return list_of_rows_to_remove


def cat_purification_cores(catNB3, catMB3, catNB6, catMB6):
    catNB3_c, catNB3_m, catNB3_r = cat_purification(catNB3)
    catMB3_c, catMB3_m, catMB3_r = cat_purification(catMB3)
    catNB6_c, catNB6_m, catNB6_r = cat_purification(catNB6)
    catMB6_c, catMB6_m, catMB6_r = cat_purification(catMB6)
    
    catNB3_HII, catNB3_HII_cand, catNB3_HII_new, catMB3_HII, catMB3_HII_cand, catMB3_HII_new, catNB6_HII, catNB6_HII_cand, catNB6_HII_new, catMB6_HII, catMB6_HII_cand, catMB6_HII_new = HII_catalogs(catNB3, catNB6, catMB3, catMB6)
    
    catNB3_HII_array = [catNB3_HII_cand,catNB3_HII,catNB3_HII_new]
    catMB3_HII_array = [catMB3_HII_cand,catMB3_HII,catMB3_HII_new]
    
    remove_list_N = remove_HII_cores(catNB3_HII_array, catNB6_m)
    remove_list_M = remove_HII_cores(catMB3_HII_array, catMB6_m)
    
    catNB6_m.remove_rows(remove_list_N)
    catMB6_m.remove_rows(remove_list_M)
    
    return catNB3_m, catMB3_m, catNB6_m,catMB6_m
    


def B6_names_to_keep_from_B3(catB6, catB3):
    keep_B6_name = [] 
    for i in range(len(catB6)):
        B3_name = catB6[i]['B3_match']
        B3row = catB3[catB3['_name'] == B3_name]
        if B3row['notes'] == 'HII?':
            if catB6[i]['notes'] == 'HII?':
                1==1
                #print('All good')
            if catB6[i]['notes'] == 'HII':
                1==1
                # Only want unconfirmed HII regions
            else:
                #print('B3 is marked HII?, while B6 is not')
                keep_B6_name += [catB6[i]['_name']]
                
    return keep_B6_name

def catalog_from_names(cat, names):
    temp = [np.nan] * len(cat)
    for i in range(len(cat)):
        if cat['_name'][i] in names:
            temp[i] = True
        else:
            temp[i] = False
    return cat[temp]

import numpy as np

def get_B6_HII_cat_from_B3(catB3_HII, catB6):
    B6_names = []
    B6_indexes = np.array([])
    for row in catB3_HII:
        B6_name = catB6[catB6['B3_match'] == row['_name']]['_name'].value
        B6_index = np.where(catB6['_name'] == B6_name)
        if B6_index[0].size > 0:
            B6_index = int(B6_index[0])
        B6_indexes = np.append(B6_indexes, B6_index)
    B6_indexes = B6_indexes.astype(int)
    return catB6[[B6_indexes]]


def HII_catalogs(catNB3, catNB6, catMB3, catMB6):
    catNB3_HII_new = catNB3[catNB3['notes'] == 'HII_new?']
    catNB3_HII_cand = catNB3[catNB3['notes'] == 'HII?']
    catNB3_HII = catNB3[catNB3['notes'] == 'HII']

    catMB3_HII_new = catMB3[catMB3['notes'] == 'HII_new?']
    catMB3_HII_cand = catMB3[catMB3['notes'] == 'HII?']
    catMB3_HII = catMB3[catMB3['notes'] == 'HII']

    
    catNB6_HII_cand = get_B6_HII_cat_from_B3(catNB3_HII_cand, catNB6)
    catNB6_HII = get_B6_HII_cat_from_B3(catNB3_HII, catNB6)
    catNB6_HII_new = get_B6_HII_cat_from_B3(catNB3_HII_new, catNB6)

    catMB6_HII_cand = get_B6_HII_cat_from_B3(catMB3_HII_cand, catMB6)
    catMB6_HII = get_B6_HII_cat_from_B3(catMB3_HII, catMB6)
    catMB6_HII_new = get_B6_HII_cat_from_B3(catMB3_HII_new, catMB6)
    
    return(catNB3_HII, catNB3_HII_cand, catNB3_HII_new, catMB3_HII, catMB3_HII_cand, catMB3_HII_new, 
           catNB6_HII, catNB6_HII_cand, catNB6_HII_new, catMB6_HII, catMB6_HII_cand, catMB6_HII_new)


def HII_catalogs_old(catNB3, catNB6, catMB3, catMB6):
    NB3names_keep_B6 = catNB6[catNB6['notes'] == 'HII?']['B3_match']
    NB3names_keep_B6 = NB3names_keep_B6[~NB3names_keep_B6.mask]
    NB3names_keep_B6 = list(set(NB3names_keep_B6))

    NB3names_keep_HIIq = catNB3[catNB3['notes'] == 'HII?']['_name']
    NB3names_keep_HIIq = list(set(NB3names_keep_HIIq))

    NB3names_keep_HII = catNB3[catNB3['notes'] == 'HII']['_name']
    NB3names_keep_HII = list(set(NB3names_keep_HII))
    
    NB3names_keep_HII_newq = catNB3[catNB3['notes'] == 'HII_new?']['_name']
    NB3names_keep_HII_newq = list(set(NB3names_keep_HII_newq))


    NB3_HII_candidates = set(NB3names_keep_B6 + NB3names_keep_HIIq)
    NB3_HII = NB3names_keep_HII
    NB3_HII_new = NB3names_keep_HII_newq


    ###########
    MB3names_keep_B6 = catMB6[catMB6['notes'] == 'HII?']['B3_match']
    MB3names_keep_B6 = MB3names_keep_B6[~MB3names_keep_B6.mask]
    MB3names_keep_B6 = list(set(MB3names_keep_B6))

    MB3names_keep_HIIq = catMB3[catMB3['notes'] == 'HII?']['_name']
    MB3names_keep_HIIq = list(set(MB3names_keep_HIIq))

    MB3names_keep_HII = catMB3[catMB3['notes'] == 'HII']['_name']
    MB3names_keep_HII = list(set(MB3names_keep_HII))

    MB3names_keep_HII_newq = catMB3[catMB3['notes'] == 'HII_new?']['_name']
    MB3names_keep_HII_newq = list(set(MB3names_keep_HII_newq))
    

    MB3_HII_candidates = set(MB3names_keep_B6 + MB3names_keep_HIIq)
    MB3_HII = MB3names_keep_HII
    MB3_HII_new = MB3names_keep_HII_newq
    
    
    
    NB6_HII_candidates = B6_names_to_keep_from_B3(catNB6, catNB3)
    MB6_HII_candidates = B6_names_to_keep_from_B3(catMB6, catMB3)

    NB6_HII = list((catNB6[catNB6['notes'] == 'HII']['_name']).data)
    MB6_HII = list((catMB6[catMB6['notes'] == 'HII']['_name']).data)
    
    
    
    
    cat_NB3_HII_candidates = catalog_from_names(catNB3, NB3_HII_candidates)
    cat_MB3_HII_candidates = catalog_from_names(catMB3, MB3_HII_candidates)
    cat_NB3_HII = catalog_from_names(catNB3, NB3_HII)
    cat_MB3_HII = catalog_from_names(catMB3, MB3_HII)
    cat_NB3_HII_new = catalog_from_names(catNB3, NB3_HII_new)
    cat_MB3_HII_new = catalog_from_names(catMB3, MB3_HII_new)
    cat_NB6_HII_candidates = catalog_from_names(catNB6, NB6_HII_candidates)
    cat_MB6_HII_candidates = catalog_from_names(catMB6, MB6_HII_candidates)
    cat_NB6_HII = catalog_from_names(catNB6, NB6_HII)
    cat_MB6_HII = catalog_from_names(catMB6, MB6_HII)
    
    return(cat_NB3_HII_candidates, cat_MB3_HII_candidates, cat_NB3_HII, cat_MB3_HII, cat_NB3_HII_new, cat_MB3_HII_new, cat_NB6_HII_candidates, 
           cat_MB6_HII_candidates, cat_NB6_HII, cat_MB6_HII)



