import pandas
def load_expr(array_exp_file, num_sel, delimiter='\t'):
"""Load the dataset in a pandas dataframe"""
    if os.path.isfile(array_exp_file):
        expr = pandas.csv_read(array_exp_file, sep=delimiter)
    (rows, cols) = expr.shape
    if num_sel > cols:
        num_sel = cols
    expr_sel = expr.ix[:,0:num_sel]
    expr_test = expr.ix[:,num_sel::]
    return expr_sel, expr_test, expr

def get_coef_var(expr, n):
"""
    Determine the coefficient of variation for each gene
    and return the top N genes
"""
    cv = expr.mean()/expr.std()
    cv_sel = cv.sort()[:n]
    return cv_sel



#Rank the genes and select the genes above the cut off
#Build a new matrix of the selected genes
#create a heatmap of that matrix
