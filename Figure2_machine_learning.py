# %%
# Usage:
# nohup python GLM_TAD.py 'ccle_ctrp_matched_final_drugsMat.txt' 'TADscore-GSVA-10kb-CCLE-CTRP_log2.csv' 'columnname_final_expression.txt' &
# nohup python GLM_TAD.py <drugMat_file> <TADScore_file> <TADScore_cols> &

# %% Imports

import pandas as pd
import numpy as np
import scipy as sp
from multiprocessing import Pool
from sklearn.model_selection import train_test_split
from sklearn.linear_model import Lasso, LassoCV
from sklearn.preprocessing import StandardScaler

import warnings
warnings.filterwarnings('ignore')



# %%

# drugMat_file = str(sys.argv[1])
# TADScore_file = str(sys.argv[2])
# TADScore_cols = str(sys.argv[3])

results = list()


print("TAD Based Analysis..")
drugMat_file = 'ccle_ctrp_matched_final_drugsMat.txt'
TADScore_file = 'TADscore-GSVA-unionTADs-CCLE-CTRP_without_white_space_4dec22.csv'
# TADScore_file = './random_TADs/TADscoreRandomboundaries-GSVA-10kb-CCLE-CTRP.csv'
TADScore_cols = 'columnname_final_expression.txt'


# print("Gene Set Basis..")
# drugMat_file = 'ccle_ctrp_matched_final_drugsMat.txt'
# TADScore_file = './Genes_set_basis/ccle_ctrp_matched_final_expression_binary_seection_qqnorm_log.csv'
# TADScore_cols = 'columnname_final_expression.txt'



# %%

def clean_df(df):
    df.replace([np.inf, -np.inf, np.nan], 0, inplace=True)
    print("dataframe got cleaned..")


def myLasso(args):
    print("input data index:",args[0])
    return [drugs[args[0]], myLasso_(args[1])]


def myLasso_(t):
    
    cols = t.columns
    # Identifying X and y
    X = t.loc[:, cols[:len(cols)-4]].astype('float64')
    y = t.loc[:, 'PIC50'].astype('float64')

    # Train Test split
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3) #, random_state=10)

    # Standardization
    scaler = StandardScaler().fit(X_train) 
    X_train = scaler.transform(X_train)
    X_test = scaler.transform(X_test)

    # Cross Validation
    model = LassoCV(cv=5, max_iter=100000, n_jobs=-1) #, random_state=0
    model.fit(X_train, y_train)

    # Best Model
    lasso_best = Lasso(alpha=model.alpha_)
    lasso_best.fit(X_train, y_train)

    # Training data
    pred_train = lasso_best.predict(X_train)
    corr_train = sp.stats.spearmanr(a=y_train, b=pred_train, axis=1)
    ct1 = round(corr_train.correlation, 2)
    cp1 = round(corr_train.pvalue, 2)

    # Test data
    pred_test = lasso_best.predict(X_test)
    corr_test = sp.stats.spearmanr(a=y_test, b=pred_test, axis=1)
    ct2 = round(corr_test.correlation, 2)
    cp2 = round(corr_test.pvalue, 2)

    return [ct1, cp1, ct1*ct1, ct2, cp2, ct2*ct2]

# %%

df_drugMat = pd.read_csv(drugMat_file, sep='\t')
clean_df(df_drugMat)


# %%

df_drugMat['PIC50'] = -np.log10((df_drugMat['IC50']+1) * (10**(-6)))
clean_df(df_drugMat)
df_drugMat.index = df_drugMat.CELLS

col_names = pd.read_csv(TADScore_cols, sep='\t', header=None).loc[:,0].values
df_TAD_Score = pd.read_csv(TADScore_file, sep='\t')
df_TAD_Score.columns = col_names
df_TAD_Score = df_TAD_Score.T

# %%

data = df_TAD_Score.join(df_drugMat)

# %%

drugs_name = data.DRUGS
drugs = list(data.DRUGS.unique())

# %%

input_data = []
for i in range(len(drugs)):
    temp = data[data.DRUGS == drugs[i]]
    input_data.append([i, temp.shape[0], temp])

 # %%


d = pd.DataFrame(input_data)
d.columns = ['i', 'pc', 'd']
d = d.sort_values(by='pc', ascending=False).iloc[400:,]
d = d.drop(columns=['pc'])
input_data = d.values.tolist()


# %%

nthreads = 16
pool = Pool(nthreads)
results = pool.map(myLasso, input_data)
print(results)
print("---------------------------")



# %%

one = []
for i, t in results:
    one.append([i] + t)
two = pd.DataFrame(one)
two.columns = ['Drugs','Corr_train', 'P_train', 'R2_train','Corr_test', 'P_test', 'R2_test']
two.to_csv('final_log_400_.csv')

# %%




