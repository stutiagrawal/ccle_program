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
    return df.sort(axis=1)

def convert_to_df_and_sort(dictionary, col_name):
    tmp = pandas.Series(dictionary)
    tmp = pandas.DataFrame(tmp, columns=[col_name])
    tmp = tmp.sort(columns=col_name)
    return tmp

def get_correlation(expr, drug, drug_name):
    rows_expr, cols_expr = expr.shape
    rows_drug, cols_drug = drug.shape
    corr = dict()
    p = dict()
    print "Num of celllines in expr: %d" %cols_expr
    print "Num of cell lines in drug: %d" %cols_drug
    if(cols_expr == cols_drug):
        for gene in expr.index:
            (pearson_coeff, p_value) = scipy.stats.pearsonr(expr.loc[gene], drug.loc[drug_name])
            corr[gene] = pearson_coeff
            p[gene] = p_value
    else:
        raise Exception("number of cell lines do not match")
    corr = convert_to_df_and_sort(corr, "pearson")
    p = convert_to_df_and_sort(p, "p_value")
    combine = pandas.concat([corr, p], axis =1)
    combine = combine.sort(columns="pearson")
    return combine

def get_updated_expr(expr, genes_sel):
    """ build a new matrix with selected genes"""

    to_be_removed = set()

    for gene_name in expr.index:
        if gene_name not in genes_sel.index:
            to_be_removed.add(gene_name)

    for gene_name in to_be_removed:
        expr = expr.drop(gene_name, axis=0)
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
    if "Accession" in expr.keys():
        expr = expr.drop("Accession", axis=1)

    drug = load_data_frame(args.drug_dataset, index_col="COMPOUND_NAME")
    #Remove cell lines which do not have drug response
    to_be_removed = set()
    for i in xrange(len(drug.loc[args.drug])):
        if drug.loc[args.drug][i] == -1:
            to_be_removed.add(drug.keys()[i])

    for cell_line in to_be_removed:
        drug = drug.drop(cell_line, axis=1)

    #cv_sel = get_coef_var(expr, args.n_features)
    print "selecting cell lines"
    drug_cell_lines = get_col_names(drug)
    expr = select_columns(drug_cell_lines, expr)
    expr_cell_lines = get_col_names(expr)
    drug = select_columns(expr_cell_lines, drug)

    print "getting correlation"
    correlation = get_correlation(expr, drug, args.drug)
    print len(correlation)
    min_corr = correlation.ix[0:args.n_features/2, :]
    max_corr = correlation.ix[len(correlation) - args.n_features/2:, :]
    sel = pandas.concat([min_corr, max_corr], axis=0)
    expr = get_updated_expr(expr, sel)
    print sel
    fp = open("cell_line.txt", "r")
    cell_line_names = dict()
    for line in fp:
        line = line.rstrip().split()
        if len(line) == 2:
            cell_line_names[line[1]] = line[0]
    expr = expr.rename(columns=cell_line_names)
    #for gene in sel.index:
    #    print "%s\t%s\t%s\n" %(gene, p.ix[gene]["p_value"], corr_sorted[gene]["pearson"])
    #labels = expr["Description"]
    #expr = pandas.DataFrame(expr, index=labels)
    expr.to_csv("heatmap.txt", sep="\t")
    #create_heatmap(expr)


