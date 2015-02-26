import pandas
import numpy
import matplotlib.pyplot as plt
import argparse
import os
import math
import scipy.stats

def load_data_frame(array_exp_file, delimiter='\t', index_col="Description"):
    """Load the dataset in a pandas dataframe"""

    if os.path.isfile(array_exp_file):
        expr = pandas.read_csv(array_exp_file, sep=delimiter)
    if "Accession" in expr.keys():
        expr = expr.drop("Accession", axis=1)

    expr = expr.drop_duplicates(index_col)
    expr.index = expr[index_col]
    expr = expr.drop(index_col, axis=1)
    return expr

def get_coef_var(expr, n):
    """
    Determine the coefficient of variation for each gene
    and return the top N genes
    """

    cv = expr.mean(axis=1)/expr.std(axis=1)
    cv.sort(ascending=False)
    print cv
    #print cv.shape, n, type(n)
    cv_sel = cv[:n]
    print cv_sel
    return cv_sel

def select_columns(col_names, df):
    """ Select specific cell lines """
    #print len(expr.keys())
    for column in df.keys():
        if column not in col_names:
            df = df.drop(column, axis=1)
    return df

def get_correlation(expr, drug, drug_name):
    rows_expr, cols_expr = expr.shape
    rows_drug, cols_drug = drug.shape
    corr = dict()
    print "Num of celllines in expr: %d" %cols_expr
    print "Num of cell lines in drug: %d" %cols_drug
    if(cols_expr == cols_drug):
        for gene in expr.index:
            (pearson_coeff, p_value) = scipy.stats.pearsonr(expr.loc[gene], drug[drug_name])
            corr[gene] = pearson_coeff
    else:
        raise Exception("number of cell lines do not match")
    return corr

def get_updated_expr(expr, cv_sel):
    """ build a new matrix with selected genes"""

    rows, cols = expr.shape
    for i in xrange(rows):
        if i not in cv_sel.keys():
            expr = expr.drop(i, axis=0)
    print expr
    return expr

#create a heatmap of that matrix

def create_heatmap(expr):
    print expr
    #fig, ax = plt.subplots()
    #plt.pcolor(numpy.random.rand(4,4))

    plt.pcolor(expr)
    #plt.yticks(np.arange(0.5, len(expr.index), 1), expr.index)
    #plt.xticks(np.arange(0.5, len(expr.columns), 1), expr.columns)
    plt.savefig("heatmap.png")
    plt.clf()

def get_col_names(df):
    all_names = set()
    for name in df.keys():
        name = name.rstrip()
        all_names.add(name)
    return all_names

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='ccle.py', description='Generate heatmaps')
    parser.add_argument("--expr", type=str, help='Path to expression dataset', required=True)
    parser.add_argument("--n_samples", type=int, help='No of samples selected. Enter 0 to select all',
                        required=True)
    parser.add_argument("--n_features", type=int, help='No of features selected', required=True)
    parser.add_argument("--drug_dataset", type=str, help='path to drug response dataset', required=True)
    parser.add_argument("--drug", type=str, help='Name of drug to correlate with', required=True)
    args = parser.parse_args()

    expr = load_data_frame(args.expr, index_col="Description")
    drug = load_data_frame(args.drug_dataset, index_col="COMPOUND_NAME")

    #Remove cell lines which do not have drug response
    for i in xrange(len(drug[args.drug])):
        if drug[args.drug][i] == -1:
            drug = drug.drop(drug.keys[i], axis=1)

    #cv_sel = get_coef_var(expr, args.n_features)
    #expr = get_updated_expr(expr, cv_sel)
    print "selecting cell lines"
    drug_cell_lines = get_col_names(drug)
    expr = select_columns(drug_cell_lines, expr)
    expr_cell_lines = get_col_names(expr)
    drug = select_columns(expr_cell_lines, drug)

    #for value in drug.keys():
    #    print value
    #print expr.ix[0:10, :]
    print "getting correlation"
    corr = get_correlation(expr, drug, args.drug)
    print corr
    #labels = expr["Description"]
    #expr = pandas.DataFrame(expr, index=labels)
    #expr.to_csv("gnf_out_1.txt", sep="\t")
    #create_heatmap(expr)


